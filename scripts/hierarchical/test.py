from tree_class import *

grid = GridParms(np.array([1, 2, 3, 4]), np.array([1, 1, 1, 1]), np.array([0.0, 0.0, 0.0, 0.0]))
myroot = Root(grid, 5, 1)
myroot.left = InternalNode(myroot, grid, 4, 2)
myroot.right = ExternalNode(myroot, grid, 3)
myroot.left.left = ExternalNode(myroot.left, grid, 4)
myroot.left.right = ExternalNode(myroot.left, grid, 5)

mytree = Tree(myroot)
# mytree.printTree()

n_string = "((4 6 7)(8 3))(5 2 9 1)"
n = np.array([4, 6, 7, 8, 3, 5, 2, 9, 1])
binsize = np.ones(9, dtype=int)
liml = np.zeros(9)
grid = GridParms(n, binsize, liml)

def buildBinaryTree(node, n_string, incr):
    sigma = 0
    if (n_string[0] == "("):
        for i, ele in enumerate(n_string):
            if (ele == "("):
                sigma += 1
            elif (ele == ")"):
                sigma -= 1

            if (sigma == 0):
                break
        
            n_string0 = n_string[1:i]

        n_string1 = n_string[i+2:-1]
        n0 = n_string[:i+1].replace("(", " ").replace(")", " ").split()
        n1 = n_string[i+1:].replace("(", " ").replace(")", " ").split()
        n0 = np.array([int(ele) for ele in n0])
        n1 = np.array([int(ele) for ele in n1])
        print(n_string0)
        print(n_string1)
        print("n0:", n0, "n1:", n1)

        grid0 = GridParms(n0, node.grid.binsize[:n0.size], node.grid.liml[:n0.size])
        grid1 = GridParms(n1, node.grid.binsize[n0.size:], node.grid.liml[n0.size:])
        node.left = InternalNode(grid0, 1, incr)
        buildBinaryTree(node.left, n_string0, incr+1)
        node.right = InternalNode(grid1, 1, incr)
        buildBinaryTree(node.right, n_string1, grid1, incr+1)
    else:
        n = n_string.split()
        n = [int(ele) for ele in n]
        print(n)
        grid = GridParms(n, node.grid.binsize[:n.size], node.grid.liml[:n.size])
        node
    return

newroot = Root(grid, 1, 0)
buildBinaryTree(newroot, n_string, grid, )