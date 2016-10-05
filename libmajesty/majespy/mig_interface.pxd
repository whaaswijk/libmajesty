from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from node_utils cimport node, nodeid

cdef extern from "xmg.h" namespace "majesty":
    cdef cppclass xmg:
        xmg() except +
        xmg(xmg & &) except +

        bool is_mig() const
        unsigned nin() const
        unsigned nnodes() const
        unsigned nout() const

        const vector[node]& nodes() const
        const vector[nodeid]& outputs() const
        const vector[bool]& outcompl() const
        bool equals(const xmg &) const

        string to_verilog() const


cdef extern from "mig_interface.h" namespace "majesty":
    cdef enum MoveType:
        MAJ3_PROP = 0,
        INVERTER_PROP,
        DIST_RIGHT_LEFT,
        SWAP_TERNARY,
        DIST_LEFT_RIGHT,
        MAJ3_XXY,
        MAJ3_XYY

    cdef struct move:
        MoveType type
        nodeid nodeid1
        nodeid nodeid2
        nodeid nodeid3

    cdef cppclass mig_manager:
        mig_manager() except +
        mig_manager(unsigned seed) except +

        unsigned get_seed()
        void set_seed(unsigned seed)

        xmg* create_random_graph(unsigned ninputs, unsigned nnodes)
        xmg* create_random_graph(unsigned ninputs, unsigned nnodes, unsigned noutputs)
        xmg* random_mig_decomposition(unsigned ninputs)

    unsigned get_nr_unary_moves()
    unsigned get_nr_binary_moves()
    unsigned get_nr_ternary_moves()
    unsigned get_nr_edge_types()

    xmg* apply_move(const xmg &, move &)

    float compute_reward(const xmg &, const xmg &)
    vector[move] compute_moves(const xmg&)

    xmg* mig_string_decompose(const string &)
    xmg* mig_expression_decompose(unsigned ninputs, const string &)
    xmg* mig_int_decompose(unsigned ninputs, unsigned truth_table)
    xmg* get_optimum_mig(const xmg&) except + # Can throw is Cirkit is not found
    xmg* get_optimum_xmg(const xmg&) except + # Can throw is Cirkit is not found
    
    xmg* strash_xmg(const xmg&)
    xmg* remove_duplicates(const xmg&)
    
    xmg* verilog_to_xmg_ptr(const string &)

cdef extern from "maj_io.h" namespace "majesty":
    xmg* ptr_read_bench(const string & filename)
    xmg* ptr_read_verilog(const string & filename)
    void write_verilog(const xmg &, const string & filename)

cdef extern from "lut_optimize.h" namespace "majesty":
    xmg* ptr_lut_area_strategy(const xmg &, unsigned, unsigned) except + # Can throw is Cirkit or Redis server are not found
    xmg* ptr_lut_area_timeout_strategy(const xmg & , unsigned, unsigned, unsigned) except + # Can throw is Cirkit or Redis server are not found
