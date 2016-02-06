
// VF_Table.h

#ifndef VF_Table_h
#define VF_Table_h

struct vf_t
{
    double v; // voltage
    double f; // frequency
};

#define VF_TABLE_SIZE 3001

const struct vf_t vf_table[VF_TABLE_SIZE];

#endif // VF_Table_h
