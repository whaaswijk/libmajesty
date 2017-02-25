
from typing import Union
cimport cython
import numpy as np
cimport numpy as np
from libcpp cimport bool as cppbool
from libcpp.string cimport string
from libcpp.vector cimport vector

from mig_interface cimport mig_manager, xmg, move, MoveType
cimport mig_interface
from node_utils cimport node, nodeid 
cimport node_utils
from graphviz import Digraph


# Temporary
cdef unsigned _nr_unary_move = mig_interface.get_nr_unary_moves()
cdef unsigned _nr_binary_move = mig_interface.get_nr_binary_moves()
cdef unsigned _nr_ternary_move = mig_interface.get_nr_ternary_moves()

_move_type_str = {mig_interface.MAJ3_PROP: 'Majority',
                  mig_interface.INVERTER_PROP: 'Inverted Propagation',
                  mig_interface.DIST_LEFT_RIGHT : 'Distributivity L->R',
                  mig_interface.SWAP_TERNARY : 'Swap3',
                  mig_interface.DIST_RIGHT_LEFT : 'Distributivity R->L',
                  mig_interface.MAJ3_XXY: 'Majority XXY',
                  mig_interface.MAJ3_XYY: 'Majority XYY',
                  mig_interface.SUBSTITUTION: 'Substitution',
                  mig_interface.RELEVANCE: 'Relevance'
                  }

_move_type_color = {mig_interface.MAJ3_PROP : 'brown1',
                    mig_interface.INVERTER_PROP : 'cadetblue1',
                    mig_interface.MAJ3_XXY : 'red',
                    mig_interface.MAJ3_XYY : 'blue',
                    mig_interface.SWAP_TERNARY : 'orange',
                    mig_interface.DIST_LEFT_RIGHT : 'green',
                    mig_interface.SUBSTITUTION : 'black',
                    mig_interface.RELEVANCE : 'yellow',
                    mig_interface.DIST_RIGHT_LEFT : 'brown'}


cdef class PyMove:
    cdef move c_move
    def __cinit__(self, move_type=None, node_id1=None, node_id2=None, node_id3=None):
        if move_type is not None:
            assert node_id1 is not None, "No node_id given"

            assert (move_type < _nr_unary_move) != ((node_id2 is not None) or (node_id3 is not None)),\
                "MoveType doesn't correspond to nb of params"
            assert (move_type < _nr_unary_move+_nr_binary_move) != (node_id3 is not None),\
                "MoveType doesn't correspond to nb of params"

            self.c_move.type = <MoveType>move_type
            self.c_move.nodeid1 = node_id1
            if node_id2 is not None:
                self.c_move.nodeid2 = node_id2
            if node_id3 is not None:
                self.c_move.nodeid3 = node_id3

    cdef set_data(self, move move_data):
        self.c_move = move_data

    def get_move_type(self):
        return self.c_move.type

    def get_involved_nodes(self):
        if self.c_move.type < _nr_unary_move:
            return [self.c_move.nodeid1]
        elif self.c_move.type < _nr_unary_move+_nr_binary_move:
            return [self.c_move.nodeid1, self.c_move.nodeid2]
        else:
            return [self.c_move.nodeid1, self.c_move.nodeid2, self.c_move.nodeid3]

    def as_tuple(self):
        if self.c_move.type < _nr_unary_move:
            return self.c_move.type, self.c_move.nodeid1
        elif self.c_move.type < _nr_unary_move+_nr_binary_move:
            return self.c_move.type, self.c_move.nodeid1, self.c_move.nodeid2
        else:
            return self.c_move.type, self.c_move.nodeid1, self.c_move.nodeid2, self.c_move.nodeid3

    def __getitem__(self, index):
        return self.as_tuple().__getitem__(index)

    def __repr__(self):
        return "Move({}:{},{})".format(self.c_move.type,
                                      _move_type_str[self.c_move.type],
                                      self.get_involved_nodes())

    def __richcmp__(self, el, int op):
        """
        See http://cython.readthedocs.io/en/latest/src/userguide/special_methods.html#rich-comparisons for explanation
        :param el: Element to be compared
        :param op: Comparison operation type
        :return: Comparison value
        """
        if op == 0:
            return self.as_tuple() < el.as_tuple()
        elif op == 1:
            return self.as_tuple() <= el.as_tuple()
        elif op == 2:
            return self.as_tuple() == el.as_tuple()
        elif op == 3:
            return self.as_tuple() != el.as_tuple()
        elif op == 4:
            return self.as_tuple() > el.as_tuple()
        elif op == 5:
            return self.as_tuple() >= el.as_tuple()
        else:
            raise NotImplementedError('op code unknown ' + str(op))

    def __hash__(self):
        return hash(self.as_tuple())

    def __getstate__(self):
        return self.c_move.type, self.c_move.nodeid1, self.c_move.nodeid2, self.c_move.nodeid3

    def __setstate__(self, data):
        self.c_move.type = <MoveType>data[0]
        self.c_move.nodeid1 = data[1]
        self.c_move.nodeid2 = data[2]
        self.c_move.nodeid3 = data[3]


cdef class PyXmg:
    cdef xmg* c_xmg
    def __cinit__(self):
        """
        Creates an empty xmg
        :return:
        """
        self.c_xmg = new xmg()

    def __dealloc__(self):
        del self.c_xmg

    cdef set_pt_to(self, xmg* pt):
        """
        Only called from cython
        :param pt: the new pointer to wrap around
        The instance is now taking care of freeing the pointer at __dealloc__ time
        """
        del self.c_xmg
        assert pt != NULL
        self.c_xmg = pt
        return self

    def get_total_nr_nodes(self) -> int:
        """
        :return: total number of nodes of the graph
        """
        return self.c_xmg.nnodes()

    def get_depth(self) -> int:
        """
        :return: the length of the critical path
        """
        return self.c_xmg.depth()    	

    def get_nr_inputs(self) -> int:
        """
        :return: number of inputs of the graph
        """
        return self.c_xmg.nin()

    def get_nr_outputs(self) -> int:
        return self.c_xmg.nout()

    def is_mig(self) -> bool:
        return self.c_xmg.is_mig()

    def get_npn_representative(self) -> PyXmg:
        cdef:
            xmg* result
        result = mig_interface.get_npn_representative(self.c_xmg[0])
        return PyXmg().set_pt_to(result)

    def get_truth_table(self) -> int:
        return mig_interface.get_truth_table(self.c_xmg[0])

    def get_outputs(self) -> ( np.ndarray, np.ndarray ):
        cdef:
            vector[nodeid] outputs
            vector[cppbool] outcompl
            np.ndarray[unsigned int, ndim=1, mode='c'] pyoutputs
            np.ndarray[unsigned int, ndim=1, mode='c'] pyoutcompl
            unsigned int num_outputs
        outputs = self.c_xmg.outputs()
        outcompl = self.c_xmg.outcompl()
        num_outputs = outputs.size()
        pyoutputs = np.empty(num_outputs, dtype=np.uint32, order='C')
        pyoutcompl = np.empty(num_outputs, dtype=np.uint32, order='C')
        for i in range(num_outputs):
            pyoutputs[i] = outputs[i]
            pyoutcompl[i] = outcompl[i]
        return pyoutputs, pyoutcompl

    def get_data(self) -> (np.ndarray, np.ndarray, np.ndarray):
        """
        Makes a copy of the underlying structure of the graph
        :return: (nodes_id, nodes_inputs, edges)
            nodes_id [N_nodes] : array of the ids of the nodes
            nodes_inputs [N_nodes, 3] : inputs ids
            edges [N_nodes, 3] : edge type
        """
        cdef:
            unsigned int num_nodes
            np.ndarray[unsigned int, ndim=1, mode='c'] nodes_id
            np.ndarray[unsigned int, ndim=2, mode='c'] nodes_inputs
            np.ndarray[int, ndim=2, mode='c'] edge_type
            vector[node] nodes
            unsigned int i
        # Ask for the data
        nodes = self.c_xmg.nodes()
        # get number of nodes of the graph
        num_nodes = nodes.size()
        # Create the numpy array
        nodes_id = np.empty(num_nodes, dtype=np.uint32, order='C')
        nodes_inputs = np.empty((num_nodes, 3), dtype=np.uint32, order='C')
        edge_type = np.empty((num_nodes, 3), dtype=np.int32, order='C')
        # Fill them
        for i in range(num_nodes):
            nodes_id[i] = nodes[i].ecrep
            nodes_inputs[i, 0] = nodes[i].in1
            nodes_inputs[i, 1] = nodes[i].in2
            nodes_inputs[i, 2] = nodes[i].in3
            if node_utils.is_pi(nodes[i]):
                edge_type[i, 0] = -1
                edge_type[i, 1] = -1
                edge_type[i, 2] = -1
            else:
                edge_type[i, 0] = node_utils.is_c1(nodes[i])
                edge_type[i, 1] = node_utils.is_c2(nodes[i])
                edge_type[i, 2] = node_utils.is_c3(nodes[i])

        return nodes_id, nodes_inputs, edge_type

    def get_adjacency_tensor(self) -> np.ndarray:
        cdef:
            int edge_type
            np.ndarray[int, ndim=3, mode='c'] adj_tensor
            vector[node] data
            node n
            node_utils.nodeid n_id
            unsigned int i
        data = self.c_xmg.nodes()
        num_nodes = self.get_total_nr_nodes()
        adj_tensor = np.zeros((2, num_nodes, num_nodes), dtype=np.int32)
        for i in range(data.size()):
            n = data[i]
            # Ignore inputs
            if node_utils.is_pi(n):
                continue
            # Update adjacency
            n_id = n.ecrep
            if node_utils.is_c1(n):
                adj_tensor[1, n_id, n.in1] += 1
            else:
                adj_tensor[0, n_id, n.in1] += 1
            if node_utils.is_c2(n):
                adj_tensor[1, n_id, n.in2] += 1
            else:
                adj_tensor[0, n_id, n.in2] += 1
            if node_utils.is_c3(n):
                adj_tensor[1, n_id, n.in3] += 1
            else:
                adj_tensor[0, n_id, n.in3] += 1
        return adj_tensor

    def get_move_inds(self):
        cdef:
            vector[move] possible_moves
            unsigned int nb_moves
            unsigned int nb_unary_moves
            unsigned int nb_binary_moves
            unsigned int nb_ternary_moves
            np.ndarray[unsigned int, ndim=2, mode='c'] unary_moves
            np.ndarray[unsigned int, ndim=2, mode='c'] binary_moves
            np.ndarray[unsigned int, ndim=2, mode='c'] ternary_moves
            unsigned int i_unary
            unsigned int i_binary
            unsigned int i_ternary
            move tmp_move

        possible_moves = mig_interface.compute_moves(self.c_xmg[0])
        nb_moves = possible_moves.size()

        nb_unary_moves = 0
        nb_binary_moves = 0
        nb_ternary_moves = 0
        for i in range(nb_moves):
            if possible_moves[i].type < _nr_unary_move:
                nb_unary_moves += 1
            elif possible_moves[i].type < _nr_unary_move+_nr_binary_move:
                nb_binary_moves += 1
            else:
                nb_ternary_moves += 1

        unary_moves = np.empty((nb_unary_moves, 2), dtype=np.uint32, order='C')
        binary_moves = np.empty((nb_binary_moves, 3), dtype=np.uint32, order='C')
        ternary_moves = np.empty((nb_ternary_moves, 4), dtype=np.uint32, order='C')
        i_unary = 0
        i_binary = 0
        i_ternary = 0
        for i in range(nb_moves):
            tmp_move = possible_moves[i]
            if tmp_move.type < _nr_unary_move:
                unary_moves[i_unary, 0] = tmp_move.type
                unary_moves[i_unary, 1] = tmp_move.nodeid1
                i_unary +=1
            elif tmp_move.type < _nr_unary_move+_nr_binary_move:
                binary_moves[i_binary, 0] = tmp_move.type
                binary_moves[i_binary, 1] = tmp_move.nodeid1
                binary_moves[i_binary, 2] = tmp_move.nodeid2
                i_binary +=1
            else:
                ternary_moves[i_ternary, 0] = tmp_move.type
                ternary_moves[i_ternary, 1] = tmp_move.nodeid1
                ternary_moves[i_ternary, 2] = tmp_move.nodeid2
                ternary_moves[i_ternary, 3] = tmp_move.nodeid3
                i_ternary +=1
        return unary_moves, binary_moves, ternary_moves

    def get_moves(self):
        cdef:
            vector[move] possible_moves
            unsigned int nb_moves

        possible_moves = mig_interface.compute_moves(self.c_xmg[0])
        nb_moves = possible_moves.size()

        move_list = list()
        for i in range(nb_moves):
            m = PyMove()
            m.set_data(possible_moves[i])
            move_list.append(m)
        return move_list

    def plot(self, highlight_action=None):
        """

        :param highlight_action: tuple of the action to be highlighted (optional)
        :return: a graphviz.dot.Digraph object.
        """
        g_ids, g_edges, g_edge_type = self.get_data()
        outputs, outcompl = self.get_outputs()
        nodes_modified = []
        action_type = None
        if highlight_action is not None:
            assert isinstance(highlight_action, PyMove), "Wrong type argument"
            nodes_modified = highlight_action.get_involved_nodes()
            action_type = highlight_action.get_move_type()
            action_color = _move_type_color[action_type]
            action_label = str(highlight_action)

        dot = Digraph()

        input_list = []
        input_subgraph = Digraph()
        input_subgraph.attr('graph', rank='same')
        output_list = []
        output_subgraph = Digraph()
        output_subgraph.attr('graph', rank='same')
        for g_id, neighbours, edge_type in zip(g_ids, g_edges, g_edge_type):
            if (edge_type == -1).all():  # Input
                input_list.append(g_id)
            else:
                dot.node(str(g_id))
                if g_id in nodes_modified:
                    dot.node(str(g_id), color=action_color, style='filled')
                for n, e_t in zip(neighbours, edge_type):
                    attrs = {}
                    if e_t==1:
                        attrs['arrowhead'] = 'dot'
                    if n in nodes_modified or g_id in nodes_modified:
                        attrs['color'] = action_color
                    dot.edge(str(n), str(g_id), **attrs)
            for i in range(len(outputs)):
                if outputs[i] == g_id:
                    attrs = {}
                    if outcompl[i] != 0:
                        attrs['arrowhead'] = 'dot'
                    dot.edge(str(g_id), 'o_'+str(i), **attrs)

        input_list = sorted(input_list)
        for g_id in input_list:
            input_subgraph.node(str(g_id), shape='doublecircle')
        # Add invisible edges to help ordering
        input_subgraph.attr('edge', style='invisible', arrowhead='none')
        for i in range(len(input_list)-1):
            input_subgraph.edge(str(input_list[i]), str(input_list[i+1]))
        dot.subgraph(input_subgraph)

        # Add invisible edges for outputs as well
        output_subgraph.attr('edge', style='invisible', arrowhead='none')
        for i, g_id in enumerate(outputs):
            output_subgraph.node('o_'+str(i), shape='doublecircle')
        for i in range(len(outputs)-1):
            output_subgraph.edge('o_'+str(i), 'o_'+str(i+1))
        dot.subgraph(output_subgraph)
        dot.graph_attr['rankdir'] = "BT"


        graph_title = "{} nodes ({}+1 inputs).".format(self.get_total_nr_nodes(), self.get_nr_inputs())
        if action_type is not None:
            graph_title += " Applying {}".format(action_label)
        dot.graph_attr['label'] = graph_title
        dot.graph_attr['labelloc'] = 't'

        return dot

    def equals(self, PyXmg py_xmg):
        return self.c_xmg.equals(py_xmg.c_xmg[0])

    def to_verilog(self):
        cdef string vstr = self.c_xmg.to_verilog()
        return vstr.decode('UTF-8')

    def __getstate__(self):
        return self.to_verilog()

    def __setstate__(self, data):
        cdef xmg* result = mig_interface.verilog_to_xmg_ptr(data.encode('UTF-8'))
        self.set_pt_to(result)

    def __hash__(self):
        return hash(str(self))

    def __str__(self):
        return self.to_verilog()

    def __richcmp__(self, el, int op):
        """
        See http://cython.readthedocs.io/en/latest/src/userguide/special_methods.html#rich-comparisons for explanation
        :param el: Element to be compared
        :param op: Comparison operation type
        :return: Comparison value
        """
        if op == 0:
            return str(self) < str(el)
        elif op == 1:
            return str(self) <= str(el)
        elif op == 2:
            return str(self) == str(el)
        elif op == 3:
            return str(self) != str(el)
        elif op == 4:
            return str(self) > str(el)
        elif op == 5:
            return str(self) >= str(el)
        else:
            raise NotImplementedError('op code unknown ' + str(op))

    def __repr__(self):
        return "PyXmg({} nodes, {})".format(self.get_total_nr_nodes(), hash(self))


cdef class MigManager:
    cdef mig_manager* c_mig_manager  # hold a C++ instance which we're wrapping
    def __cinit__(self, seed=None):
        if seed is None:
            self.c_mig_manager = new mig_manager()
        else:
            self.c_mig_manager = new mig_manager(seed)

    def __dealloc__(self):
        del self.c_mig_manager

    def get_seed(self) -> int:
        return self.c_mig_manager.get_seed()
    def set_seed(self, unsigned int s):
        self.c_mig_manager.set_seed(s)

    def create_random_graph(self, unsigned int num_input, unsigned int num_nodes) -> PyXmg:
        assert num_input + 1 <= num_nodes, "num_input can not be larger that num_nodes"
        return PyXmg().set_pt_to(self.c_mig_manager.create_random_graph(num_input, num_nodes))

    def random_mig_decomposition(self, unsigned int ninputs, strash=True) -> PyXmg:
        cdef:
            xmg* result
            xmg* strashed_result
        result = self.c_mig_manager.random_mig_decomposition(ninputs)
        if strash:
            strashed_result = mig_interface.strash_xmg(result[0])
            strashed_result = mig_interface.remove_duplicates(strashed_result[0])
            del result
            return PyXmg().set_pt_to(strashed_result)
        else:
            return PyXmg().set_pt_to(result)


def get_nr_unary_moves() -> int:
    return mig_interface.get_nr_unary_moves()
def get_nr_binary_moves() -> int:
    return mig_interface.get_nr_binary_moves()
def get_nr_ternary_moves() -> int:
    return mig_interface.get_nr_ternary_moves()
def get_nr_edge_types() -> int:
    return mig_interface.get_nr_edge_types()

def apply_move(PyXmg py_xmg, PyMove move, remove_duplicates=True, strash=True, strash_if_substitution=True) -> Union[None, PyXmg]:
    cdef:
        xmg* result
        xmg* new_result
    for node_id in move.get_involved_nodes():
        if node_id not in range(py_xmg.get_total_nr_nodes()):
            return None
    result = mig_interface.apply_move(py_xmg.c_xmg[0], move.c_move)
    if result == NULL:
        return None
    else:
        if strash or (strash_if_substitution and move.get_move_type() == mig_interface.SUBSTITUTION):
            new_result = mig_interface.strash_xmg(result[0])
            del result
            return PyXmg().set_pt_to(new_result)
        elif remove_duplicates:
            new_result = mig_interface.remove_duplicates(result[0])
            del result
            return PyXmg().set_pt_to(new_result)
        else:
            return PyXmg().set_pt_to(result)

def mig_string_decompose(tt) -> PyXmg:
    cdef:
        xmg* result
        string truth_table
    truth_table = tt.encode('UTF-8')
    result = mig_interface.mig_string_decompose(truth_table)
    return PyXmg().set_pt_to(result)

def mig_expression_decompose(py_ninputs, py_expr) -> PyXmg:
    cdef:
        xmg* result
        unsigned int ninputs
        string expr
    ninputs = py_ninputs
    expr = py_expr.encode('UTF-8')
    result = mig_interface.mig_expression_decompose(ninputs, expr)
    return PyXmg().set_pt_to(result)

def mig_int_decompose(nin, tt) -> PyXmg:
    cdef:
        xmg* result
        unsigned int truth_table
        unsigned int ninputs
    truth_table = tt
    ninputs = nin
    result = mig_interface.mig_int_decompose(ninputs, truth_table)
    return PyXmg().set_pt_to(result)

def get_optimum_mig(PyXmg py_xmg) -> PyXmg:
    cdef:
        xmg* result
    result = mig_interface.get_optimum_mig(py_xmg.c_xmg[0])
    return PyXmg().set_pt_to(result)

def get_depth_optimum_mig(PyXmg py_xmg) -> PyXmg:
    cdef:
        xmg* result
    result = mig_interface.get_depth_optimum_mig(py_xmg.c_xmg[0])
    return PyXmg().set_pt_to(result)

def get_optimum_xmg(PyXmg py_xmg) -> PyXmg:
    cdef:
        xmg* result
    result = mig_interface.get_optimum_xmg(py_xmg.c_xmg[0])
    return PyXmg().set_pt_to(result)

def get_npn_representative(PyXmg py_xmg) -> PyXmg:
    cdef:
        xmg* result
    result = mig_interface.get_npn_representative(py_xmg.c_xmg[0])
    return PyXmg().set_pt_to(result)

def resyn2(PyXmg py_xmg) -> PyXmg:
    cdef:
        xmg* result
    result = mig_interface.resyn2(py_xmg.c_xmg[0])
    return PyXmg().set_pt_to(result)

def get_npn_representative(PyXmg py_xmg) -> int:
    return mig_interface.get_truth_table(py_xmg.c_xmg[0])

def strash_xmg(PyXmg py_xmg) -> PyXmg:
    cdef:
        xmg* result
    result = mig_interface.strash_xmg(py_xmg.c_xmg[0])
    return PyXmg().set_pt_to(result)

def remove_duplicates(PyXmg py_xmg) -> PyXmg:
    cdef:
        xmg* result
    result = mig_interface.remove_duplicates(py_xmg.c_xmg[0])
    return PyXmg().set_pt_to(result)

def compute_reward(PyXmg xmg_initial, PyXmg xmg_final) -> float:
    return mig_interface.compute_reward(xmg_initial.c_xmg[0], xmg_final.c_xmg[0])

def read_bench(py_filename) -> PyXmg:
    cdef:
        string filename
        xmg* result
    filename = py_filename.encode('UTF-8')
    result = mig_interface.ptr_read_bench(filename)
    return PyXmg().set_pt_to(result)

def read_verilog(py_filename) -> PyXmg:
    cdef:
        string filename
        xmg* result
    filename = py_filename.encode('UTF-8')
    result = mig_interface.ptr_read_verilog(filename)
    return PyXmg().set_pt_to(result)

def write_verilog(PyXmg py_xmg, py_filename) -> None:
    cdef:
        string filename
    filename = py_filename.encode('UTF-8')
    mig_interface.write_verilog(py_xmg.c_xmg[0], filename)

def verilog_to_xmg(str):
    cdef xmg* result = mig_interface.verilog_to_xmg_ptr(str.encode('UTF-8'))
    return PyXmg().set_pt_to(result)

