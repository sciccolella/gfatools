//#define GFA_EXT_SEG
#ifndef __GFA_EXT_H__
#define __GFA_EXT_H__
#include <stdint.h>
#include <stdio.h>

/* Delete this define */
#define GFA_EXT_ALL

#ifdef GFA_EXT_ALL
#define GFA_EXT_L
#define GFA_EXT_S
#define GFA_EXT_W
#endif

#if (defined(GFA_EXT_S) || defined(GFA_EXT_L)) && !defined(GFA_EXT_AUX)
#define GFA_EXT_AUX
#endif

#ifdef GFA_EXT_AUX
typedef struct {
	uint32_t m_aux, l_aux;
	uint8_t *aux;
} gfa_aux_t;

typedef struct {
	uint32_t start, end; // start: starting vertex in the string graph; end: ending vertex
	uint32_t len_comp, dummy; // len_comp: the length of the complement unitig
	uint32_t m, n; // number of reads
	uint64_t *a; // list of reads
	uint64_t *r; // start and end on each read
	char **name;
} gfa_utg_t;
#endif

#ifdef GFA_EXT_L
typedef struct {
    uint64_t v_lv; // higher 32 bits: vertex_id; lower 32 bits: lv; packed together for sorting
    uint32_t w;
    int32_t rank;
    int32_t ov, ow;
    char *exons;
    uint64_t link_id:61, strong:1, del:1, comp:1; // link_id: a pair of dual arcs are supposed to have the same link_id
} gfa_arc_t;
#endif

#ifdef GFA_EXT_S
typedef struct {
    int32_t len;
    uint32_t del:16, circ:16;
    int32_t snid; // stable name ID
    int32_t soff; // stable start position
    int32_t rank; // stable rank
    int32_t cov;
    char *name, *seq;
    gfa_utg_t *utg;
    gfa_aux_t aux;
} gfa_seg_t;
#endif

#ifdef GFA_EXT_W
typedef struct {
	const char *sample;
	int32_t snid;
	int32_t hap, n_v;
	int64_t st, en;
	uint32_t *v;
	gfa_aux_t aux;
} gfa_walk_t;
#endif

// #include "gfa.h"
// #include "gfatools/gfa-priv.h"

int gfa_aux_parse(char *s, uint8_t **data, int *max);

#endif
