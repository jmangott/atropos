#ifndef TREE_CLASS_HPP
#define TREE_CLASS_HPP

#include <iostream>
#include <vector>

#include <generic/matrix.hpp>
#include <generic/storage.hpp>
#include <generic/tensor.hpp>
#include <generic/timer.hpp>
#include <generic/tree.hpp>
#include <lr/coefficients.hpp>
#include <lr/lr.hpp>

#include "coeff_class.hpp"
#include "grid_class.hpp"
#include "matrix.hpp"

struct cme_node : virtual Ensign::node<double> {
    std::array<cme_node*, 2> child;
    const grid_parms grid;
    cme_coeff coefficients;

    cme_node(const std::string _id, cme_node* const _parent, const grid_parms _grid,
             const Index _r_in, const Index _n_basisfunctions)
        : Ensign::node<double>(_id, _parent, {nullptr, nullptr}, _r_in,
                               _n_basisfunctions),
          child({nullptr, nullptr}), grid(_grid), coefficients(_grid.n_reactions, _r_in)
    {
    }
    void CalculateAB_bar(const Ensign::Matrix::blas_ops& blas);
};

struct cme_internal_node : cme_node, Ensign::internal_node<double> {
    cme_internal_node(const std::string _id, cme_internal_node* const _parent,
                      const grid_parms _grid, const Index _r_in,
                      const std::array<Index, 2> _r_out, const Index _n_basisfunctions)
        : Ensign::node<double>(_id, _parent, {nullptr, nullptr}, _r_in,
                               _n_basisfunctions),
          cme_node(_id, _parent, _grid, _r_in, _n_basisfunctions),
          Ensign::internal_node<double>(_id, _parent, _r_in, _r_out, _n_basisfunctions)
    {
    }
    void Initialize(int ncid) override;

    template <Index id> void CalculateAB(const Ensign::Matrix::blas_ops& blas);
};

struct cme_external_node : cme_node, Ensign::external_node<double> {
    std::vector<std::vector<double>> propensity;

    cme_external_node(const std::string _id, cme_internal_node* const _parent,
                      const grid_parms _grid, const Index _r_in,
                      const Index _n_basisfunctions)
        : Ensign::node<double>(_id, _parent, {nullptr, nullptr}, _r_in,
                               _n_basisfunctions),
          cme_node(_id, _parent, _grid, _r_in, _n_basisfunctions),
          Ensign::external_node<double>(_id, _parent, _grid.dx, _r_in,
                                        _n_basisfunctions),
          propensity(_grid.n_reactions)
    {
    }
    void Initialize(int ncid) override;
};

Ensign::multi_array<double, 2> CalculateKDot(const Ensign::multi_array<double, 2>& K,
                                             const cme_external_node* const node,
                                             const Ensign::Matrix::blas_ops& blas);

Ensign::multi_array<double, 2> CalculateSDot(const Ensign::multi_array<double, 2>& S,
                                             const cme_node* const node,
                                             const Ensign::Matrix::blas_ops& blas);

Ensign::multi_array<double, 2> CalculateQDot(const Ensign::multi_array<double, 2>& Qmat,
                                             const cme_internal_node* const node,
                                             const Ensign::Matrix::blas_ops& blas);

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
                             const Ensign::Matrix::blas_ops& blas) const;
    void InitializeAB_barHelper(cme_node* const node,
                                const Ensign::Matrix::blas_ops& blas) const;
    std::vector<double> NormalizeHelper(cme_node const* const node) const;

  public:
    void Read(const std::string fn);
    void Write(const std::string, const double t, const double tau,
               const double dm) const;
    void Orthogonalize(const Ensign::Matrix::blas_ops& blas) const;
    void InitializeAB_bar(const Ensign::Matrix::blas_ops& blas) const;
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

template <Index id>
void cme_internal_node::CalculateAB(const Ensign::Matrix::blas_ops& blas)
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
