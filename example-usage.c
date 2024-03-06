#include "gfa.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv) {
  gfa_t *ingfa;

  ingfa = gfa_read(argv[1]);

  printf("------ segments -----\n");
  // Iterate over the segments
  for (int i = 0; i < ingfa->n_seg; ++i) {
    printf("name: %s\t", ingfa->seg[i].name);
    printf("len: %8d\t", ingfa->seg[i].len);
    printf("seq[:5]: %.5s\t", ingfa->seg[i].seq);
    printf("n_seg[i]: %5d\t", i);
    printf("name2id: %5d\t", gfa_name2id(ingfa, ingfa->seg[i].name));
    printf("rank: %4d\t", ingfa->seg[i].rank);
    printf("[ext]cov: %3d\t", ingfa->seg[i].cov);
    printf("\n");
  }

  printf("------ vertices -----\n");
  uint32_t nvtx = gfa_n_vtx(ingfa);
  printf("gfa_n_vtx(): %d\n", nvtx);

  for (int vix = 0; vix < nvtx; ++vix) {
    int nv = gfa_arc_n(ingfa, vix);
    gfa_arc_t *av = gfa_arc_a(ingfa, vix);

    if (nv <= 0)
      continue;
    printf("vix: %d\t nv: %d\n", vix, nv);

    for (int i = 0; i < nv; ++i) {
      gfa_arc_t *avi = &av[i];
      printf("\t- h:%d [%d, %s] -> t:%d [%d, %s]\n", gfa_arc_head(*avi),
             gfa_arc_head(*avi) >> 1,
             ingfa->seg[gfa_arc_head(*avi) >> 1]
                     .name, // index is also: avi->v_lv>>33
             gfa_arc_tail(*avi), gfa_arc_tail(*avi) >> 1,
             ingfa->seg[gfa_arc_tail(*avi) >> 1].name // index is also avi->w>>1
      );
    }
  }

  printf("------ arcs -----\n");
  char *segname = argc > 2 ? argv[2] : "s1";

  // get segment from name:
  int32_t segid = gfa_name2id(ingfa, segname);
  if (segid == -1) {
    fprintf(stderr, "segment '%s' not found.\n", segname);
    return EXIT_FAILURE;
  }

  /* get the arcs of the segment
   * the segment `segid` is associated to two vertices:
   * 1. segid +
   * v+ = segid << 1 (& 0x0)
   * 2. segid-
   * v- = segid << 1 & 1
   */
  uint32_t vid = segid << 1;

  // tags or the arc
  const gfa_aux_t *arc_aux;

  /* get the outgoing arcs of `segid`
   * segid+ -> vid
   * `gfa_arc_a(GFA, vid)` is the array of arcs
   * -- (of length `gfa_arc_n(GFA, vid)`)  --
   * that leaves from vid
   */
  printf("[Segment %s (id:%d)] Outgoing:\n", segname, segid);
  gfa_arc_t *outa_vid = gfa_arc_a(ingfa, vid);
  uint32_t n_outa_vid = gfa_arc_n(ingfa, vid);
  for (int j = 0; j < n_outa_vid; ++j) {

    printf("\t%s -> %s (%c)\t", ingfa->seg[outa_vid[j].v_lv >> 33].name,
           ingfa->seg[outa_vid[j].w >> 1].name,
           "+-"[outa_vid[j].v_lv >> 32 & 1]);
    printf("link_id:%lu\t", (uint64_t) outa_vid[j].link_id);
    
#ifdef GFA_EXT_L
    printf("[ext]ex:%s\t", outa_vid[j].genes);
#endif
    /* Alternatively it is possible to write
     *
     * printf("\t%s -> %s\n", ingfa->seg[gfa_arc_head(outa_vid[j]) >> 1].name,
     *        ingfa->seg[gfa_arc_tail(outa_vid[j]) >> 1].name);
     */



    /* Get tags, this kinda complicated.
     * So, basically tags are saved as follows; suppose we have the following
     * list of `tags`
     * RC:i:40  AN:Z:annotation CV:i:33 TR:Z:transcript
     *
     * - if the tag is an `i` it is saved as
     *        `RCi((int32_t)40)`
     *   meaning that we have the tag, its format accounting
     *   for 3 uint8_t (char) + 4 bytes representing the integer
     *
     * - if the tag is a `Z` then it is saved as
     *        `ANZannotation\0`
     *   so the columns are omitted and the string as is.
     *
     *
     * Then we concatenate all the tags
     * resulting in:
     * RCi[40=0x00000028]ANZannotation\0CVi[33=0x00000021]TRZtranscript\0
     *
     *
     * Now retrieving such information is not very efficient because we need
     * to iterate over the "string" in `gfa_aux_t->aux` until we find the tag
     * we need. If such information is needed extensively,
     * it is better to write a specific data structure to access
     * directly the information needed using the same `link_id` as index.
     *
     * Note that this is what is already done in the code
     * where, for example, the tags `SR`, `L1` and `L2`
     * as saved directly in `arc->rank`, `seg[v>>1].len`, `seg[w>>1].len`
     * and the deleted from `gfa_aux_t->aux`.
     * See https://github.com/lh3/gfatools/blob/c31be8a62efc6bdea4576029f7fbe84f345a6eed/gfa-io.c#L241
     */


    arc_aux = outa_vid[j].link_id < ingfa->n_arc ? &ingfa->link_aux[outa_vid[j].link_id] : 0;
    if (arc_aux && arc_aux->l_aux) {

      uint8_t *tag_i = gfa_aux_get(arc_aux->l_aux, arc_aux->aux, "RC");
      if (tag_i) printf("%d\t", *(int32_t*)(tag_i + 1));

      uint8_t *tag_Z = gfa_aux_get(arc_aux->l_aux, arc_aux->aux, "EX");
      if (tag_Z) printf("%s\t", tag_Z + 1);

    }
    printf("\n");
  }

  /* get the incoming arcs of `segid`
   * segid- -> vid & 1 = vid++
   * `gfa_arc_a(GFA, vid)` is the array of arcs
   * -- (of length `gfa_arc_n(GFA, vid)`)  --
   * that leaves from vid
   */
  printf("[Segment %s (id:%d)] Incoming:\n", segname, segid);
  ++vid;
  gfa_arc_t *inca_vid = gfa_arc_a(ingfa, vid);
  uint32_t n_inca_vid = gfa_arc_n(ingfa, vid);
  for (int j = 0; j < n_inca_vid; ++j) {
    printf("\t%s <- %s (%c)\t", ingfa->seg[inca_vid[j].v_lv >> 33].name,
           ingfa->seg[inca_vid[j].w >> 1].name,
           "+-"[inca_vid[j].v_lv >> 32 & 1]);
    printf("link_id:%lu\t", (uint64_t) inca_vid[j].link_id);
    printf("\n");
  }

  gfa_destroy(ingfa);

  return EXIT_SUCCESS;
}
