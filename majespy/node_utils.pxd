from libc cimport stdint
from libcpp cimport bool

cdef extern from "xmg.h":
    ctypedef stdint.uint32_t nodeid
    ctypedef bool cppbool

cdef extern from "xmg.h" namespace "majesty":
    cdef struct node:
        nodeid in1, in2, in3
        # Bit 0: is primary input
        # Bit 1: is primary output
        # Bit 2: in1 is complemented
        # Bit 3: in2 is complemented
        # Bit 4: in3 is complemented
        # Bit 5: is complement of eq. class representative
        # Bit 6: flag for misc purposes
        # Bit 7: if set this is an XOR node
        stdint.uint8_t flag
        stdint.int32_t ecrep, ecnext

    bool is_pi(const node& n)
    bool is_po(const node& n)
    bool is_c1(const node& n)
    bool is_c2(const node& n)
    bool is_c3(const node& n)
    bool is_c(const node& n)
    bool is_flag_set(const node& n)
    bool is_xor(const node& n)
    bool is_and(const node& n)
    bool is_or(const node& n)
    bool is_maj(const node& n)

    cdef struct edge:
        nodeid i, j
        bool is_complemented
        bool is_virtual
