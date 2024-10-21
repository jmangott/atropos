#ifndef TREE_CLASS_HPP
#define TREE_CLASS_HPP

#include <iostream>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <generic/timer.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include <netcdf.h>

#include "coeff_class.hpp"
#include "grid_class.hpp"
#include "matrix.hpp"
#include "netcdf_check.hpp"

#ifdef __OPENMP__
#pragma omp declare reduction(+ : Ensign::multi_array<double, 2> : omp_out += omp_in)  \
    initializer(omp_priv = decltype(omp_orig)(omp_orig))
#endif

#ifdef __OPENMP__
#pragma omp declare reduction(+ : Ensign::multi_array<double, 4> : omp_out += omp_in)  \
    initializer(omp_priv = decltype(omp_orig)(omp_orig))
#endif

// General classes for the hierarchical DLR approximation
// TODO: introduce a template parameter `N` for arbitrary many outgoing legs
template <class T> struct node {
    const std::string id;

    node* const parent;
    std::array<node*, 2> child;
    const Index n_basisfunctions;

    Ensign::multi_array<T, 2> S;

    node() = default;

    node(const std::string _id, node* const _parent, std::array<node*, 2> _child,
         const Index _r_in, const Index _n_basisfunctions)
        : id(_id), parent(_parent), child(_child), n_basisfunctions(_n_basisfunctions),
          S({_r_in, _r_in})
    {
        assert(n_basisfunctions <= _r_in);
    }

    virtual ~node() = default;

    virtual bool IsInternal() const = 0;
    virtual bool IsExternal() const = 0;
    virtual void Initialize(int ncid) = 0;

    Index RankIn() const { return S.shape()[0]; }
};

template <class T> struct internal_node : virtual node<T> {
    Ensign::multi_array<T, 3> Q;
    Ensign::multi_array<T, 3> G;

    internal_node(const std::string _id, internal_node* const _parent,
                  const Index _r_in, const std::array<Index, 2> _r_out,
                  const Index _n_basisfunctions)
        : node<T>(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions),
          Q({_r_out[0], _r_out[1], _r_in}), G({_r_out[0], _r_out[1], _r_in})
    {
    }

    bool IsInternal() const override { return true; }

    bool IsExternal() const override { return false; }

    std::array<Index, 2> RankOut() const
    {
        return array<Index, 2>({Q.shape()[0], Q.shape()[1]});
    }

    void Initialize(int ncid) override;

    void Write(int ncid, int id_r_in, std::array<int, 2> id_r_out) const;

    Ensign::multi_array<T, 2> Orthogonalize(const T weight,
                                            const Ensign::blas_ops& blas);
};

template <class T> struct external_node : virtual node<T> {
    Ensign::multi_array<T, 2> X;

    external_node(const std::string _id, internal_node<T>* const _parent,
                  const Index _dx, const Index _r_in, const Index _n_basisfunctions)
        : node<T>(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions),
          X({_dx, _r_in})
    {
    }

    bool IsInternal() const override { return false; }

    bool IsExternal() const override { return true; }

    Index ProblemSize() const { return X.shape()[0]; }

    void Initialize(int ncid) override;

    void Write(int ncid, int id_r_in, int id_dx) const;

    Ensign::multi_array<T, 2> Orthogonalize(const T weight,
                                            const Ensign::blas_ops& blas);
};

struct cme_node : virtual node<double> {
    std::array<cme_node*, 2> child;
    const grid_parms grid;
    cme_coeff coefficients;

    cme_node(const std::string _id, cme_node* const _parent, const grid_parms _grid,
             const Index _r_in, const Index _n_basisfunctions)
        : node<double>(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions),
          child({nullptr, nullptr}), grid(_grid), coefficients(_grid.n_reactions, _r_in)
    {
    }
    void CalculateAB_bar(const Ensign::blas_ops& blas);
};

struct cme_internal_node : cme_node, internal_node<double> {
    cme_internal_node(const std::string _id, cme_internal_node* const _parent,
                      const grid_parms _grid, const Index _r_in,
                      const std::array<Index, 2> _r_out, const Index _n_basisfunctions)
        : node(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions),
          cme_node(_id, _parent, _grid, _r_in, _n_basisfunctions),
          internal_node<double>(_id, _parent, _r_in, _r_out, _n_basisfunctions)
    {
    }
    void Initialize(int ncid) override;

    template <Index id> void CalculateAB(const Ensign::blas_ops& blas);
};

struct cme_external_node : cme_node, external_node<double> {
    std::vector<std::vector<double>> propensity;

    cme_external_node(const std::string _id, cme_internal_node* const _parent,
                      const grid_parms _grid, const Index _r_in,
                      const Index _n_basisfunctions)
        : node(_id, _parent, {nullptr, nullptr}, _r_in, _n_basisfunctions),
          cme_node(_id, _parent, _grid, _r_in, _n_basisfunctions),
          external_node<double>(_id, _parent, _grid.dx, _r_in, _n_basisfunctions),
          propensity(_grid.n_reactions)
    {
    }
    void Initialize(int ncid) override;
};

Ensign::multi_array<double, 2> CalculateKDot(const Ensign::multi_array<double, 2>& K,
                                             const cme_external_node* const node,
                                             const Ensign::blas_ops& blas);

Ensign::multi_array<double, 2> CalculateSDot(const Ensign::multi_array<double, 2>& S,
                                             const cme_node* const node,
                                             const Ensign::blas_ops& blas);

Ensign::multi_array<double, 2> CalculateQDot(const Ensign::multi_array<double, 2>& Qmat,
                                             const cme_internal_node* const node,
                                             const Ensign::blas_ops& blas);

struct cme_lr_tree {
    cme_internal_node* root;
    std::string partition_str;
    std::vector<std::string> species_names;

    friend std::ostream& operator<<(std::ostream& os, cme_lr_tree const& tree)
    {
        tree.PrintHelper(os, tree.root);
        return os;
    }

  private:
    void PrintHelper(std::ostream& os, cme_node const* const node) const;
    void OrthogonalizeHelper(cme_internal_node* const node,
                             const Ensign::blas_ops& blas) const;
    void InitializeAB_barHelper(cme_node* const node,
                                const Ensign::blas_ops& blas) const;
    std::vector<double> NormalizeHelper(cme_node const* const node) const;

  public:
    void Read(const std::string fn);
    void Write(const std::string, const double t, const double tau,
               const double dm) const;
    void Orthogonalize(const Ensign::blas_ops& blas) const;
    void InitializeAB_bar(const Ensign::blas_ops& blas) const;
    double Normalize() const;
};

namespace WriteHelpers {
void WritePartitionStr(int ncid, const std::string partition_str);
void WriteSpeciesNames(int ncid, const std::vector<std::string> species_names);
void WriteGridParms(int ncid, const grid_parms grid);
void WriteNode(int ncid, cme_node const* const node);
} // namespace WriteHelpers

namespace ReadHelpers {
std::string ReadPartitionStr(int ncid);
std::vector<std::string> ReadSpeciesNames(int ncid);
grid_parms ReadGridParms(int ncid);
std::array<Index, 2> ReadRankOut(int ncid);
Index ReadNBasisfunctions(int ncid);
std::vector<std::vector<double>> ReadPropensity(int ncid, const Index n_reactions);
cme_node* ReadNode(int ncid, const std::string id, cme_internal_node* const parent_node,
                   const Index r_in);
} // namespace ReadHelpers

template <class T>
Ensign::multi_array<T, 2> internal_node<T>::Orthogonalize(const T weight,
                                                          const Ensign::blas_ops& blas)
{
    Ensign::multi_array<T, 2> Qmat({Ensign::prod(RankOut()), node<T>::RankIn()});
    Ensign::multi_array<T, 2> Q_R({node<T>::RankIn(), node<T>::RankIn()});
    Ensign::Tensor::matricize<2>(Q, Qmat);
    Q_R = Ensign::Tensor::ortho(Qmat, node<T>::n_basisfunctions, weight, blas);
    Ensign::Tensor::tensorize<2>(Qmat, Q);

    return Q_R;
};

template <class T>
Ensign::multi_array<T, 2> external_node<T>::Orthogonalize(const T weight,
                                                          const Ensign::blas_ops& blas)
{
    Ensign::multi_array<T, 2> X_R({node<T>::RankIn(), node<T>::RankIn()});
    X_R = Ensign::Tensor::ortho(X, node<T>::n_basisfunctions, weight, blas);

    return X_R;
};

template <Index id> void cme_internal_node::CalculateAB(const Ensign::blas_ops& blas)
{
    const Index id_c = (id == 0) ? 1 : 0;
    Index rank_out = RankOut()[id];
    Index rank_out_c = RankOut()[id_c];

#ifdef __OPENMP__
#pragma omp parallel
#endif
    {
        Ensign::multi_array<double, 2> GA_tilde({rank_out_c * RankIn(), rank_out});
        Ensign::multi_array<double, 2> Ga_tilde({rank_out_c * RankIn(), rank_out});

        Ensign::multi_array<double, 2> GA_mat_temp({rank_out * RankIn(), rank_out_c});
        Ensign::multi_array<double, 2> Ga_mat_temp({rank_out * rank_out_c, RankIn()});
        Ensign::multi_array<double, 2> GA_mat(GA_mat_temp.shape());
        Ensign::multi_array<double, 2> Ga_mat(Ga_mat_temp.shape());

        Ensign::multi_array<double, 3> GA(G.shape());
        Ensign::multi_array<double, 3> Ga(G.shape());

        Ensign::Tensor::matricize<id_c>(G, GA_mat_temp);
        Ensign::Tensor::matricize<2>(G, Ga_mat_temp);

#ifdef __OPENMP__
#pragma omp for
#endif
        for (Index mu = 0; mu < grid.n_reactions; ++mu) {
            blas.matmul(GA_mat_temp, child[id_c]->coefficients.A_bar[mu], GA_mat);
            Ensign::Tensor::tensorize<id_c>(GA_mat, GA);
            Ensign::Tensor::matricize<id>(GA, GA_tilde);
            blas.matmul_transb(Ga_mat_temp, coefficients.A[mu], Ga_mat);
            Ensign::Tensor::tensorize<2>(Ga_mat, Ga);
            Ensign::Tensor::matricize<id>(Ga, Ga_tilde);
            blas.matmul_transa(GA_tilde, Ga_tilde, child[id]->coefficients.A[mu]);

            blas.matmul(GA_mat_temp, child[id_c]->coefficients.B_bar[mu], GA_mat);
            Ensign::Tensor::tensorize<id_c>(GA_mat, GA);
            Ensign::Tensor::matricize<id>(GA, GA_tilde);
            blas.matmul_transb(Ga_mat_temp, coefficients.B[mu], Ga_mat);
            Ensign::Tensor::tensorize<2>(Ga_mat, Ga);
            Ensign::Tensor::matricize<id>(Ga, Ga_tilde);
            blas.matmul_transa(GA_tilde, Ga_tilde, child[id]->coefficients.B[mu]);
        }
    }
};

#endif
