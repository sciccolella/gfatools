#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "gfa.h"
#include "gfa-priv.h"

int gfa_aux_parse(char *s, uint8_t **data, int *max);

int gfa_parse_S(gfa_t *g, char *s) {
  int i, is_ok = 0;
  char *p, *q, *seg = 0, *seq = 0, *rest = 0;
  uint32_t sid, len = 0;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) seg = q;
      else if (i == 1) {
        seq = q[0] == '*' ? 0 : gfa_strdup(q);
        is_ok = 1, rest = c ? p + 1 : 0;
        break;
      }
      ++i, q = p + 1;
      if (c == 0) break;
    }
  }
  if (is_ok) { // all mandatory fields read
    int l_aux, m_aux = 0, LN = -1;
    uint8_t *aux = 0, *s_LN = 0;
    gfa_seg_t *s;
    l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
    s_LN = l_aux ? gfa_aux_get(l_aux, aux, "LN") : 0;
    if (s_LN && s_LN[0] == 'i') {
      LN = *(int32_t *) (s_LN + 1);
      l_aux = gfa_aux_del(l_aux, aux, s_LN);
    }
    if (seq == 0) {
      if (LN >= 0) len = LN;
    } else len = strlen(seq);
    if (LN >= 0 && len != LN && gfa_verbose >= 2)
      fprintf(stderr, "[W] for segment '%s', LN:i:%d tag is different from sequence length %d\n", seg, LN, len);
    sid = gfa_add_seg(g, seg);
    s = &g->seg[sid];
    s->len = len, s->seq = seq;
    if (l_aux > 0) {
      uint8_t *s_SN = 0, *s_SO = 0, *s_SR = 0;
      s_SN = gfa_aux_get(l_aux, aux, "SN");
      if (s_SN && *s_SN == 'Z') { // then parse stable tags
        s->snid = gfa_sseq_add(g, (char *) (s_SN + 1)), s->soff = 0;
        l_aux = gfa_aux_del(l_aux, aux, s_SN);
        s_SO = gfa_aux_get(l_aux, aux, "SO");
        if (s_SO && *s_SO == 'i') {
          s->soff = *(int32_t *) (s_SO + 1);
          l_aux = gfa_aux_del(l_aux, aux, s_SO);
        }
      }
      s_SR = gfa_aux_get(l_aux, aux, "SR");
      if (s_SR && *s_SR == 'i') {
        s->rank = *(int32_t *) (s_SR + 1);
        if (s->rank > g->max_rank) g->max_rank = s->rank;
        l_aux = gfa_aux_del(l_aux, aux, s_SR);
      }
#ifdef GFA_EXT_S
      uint8_t *s_CV = gfa_aux_get(l_aux, aux, "CV");
      if (s_CV && *s_CV == 'i') {
        s->cov = *(int32_t*)(s_CV + 1);
        l_aux = gfa_aux_del(l_aux, aux, s_CV);
      } else {
        s->cov = 0;
      }
#endif
      gfa_sseq_update(g, s);
    }
    if (l_aux > 0)
      s->aux.m_aux = m_aux, s->aux.l_aux = l_aux, s->aux.aux = aux;
    else if (aux)
      free(aux);
  } else return -1;
  return 0;
}

int gfa_parse_L(gfa_t *g, char *s) {
//  fprintf(stderr, "[%s] gfa-ext\n", __func__);
  int i, oriv = -1, oriw = -1, is_ok = 0;
  char *p, *q, *segv = 0, *segw = 0, *rest = 0;
  int32_t ov = INT32_MAX, ow = INT32_MAX;
  for (i = 0, p = q = s + 2;; ++p) {
    if (*p == 0 || *p == '\t') {
      int c = *p;
      *p = 0;
      if (i == 0) {
        segv = q;
      } else if (i == 1) {
        if (*q != '+' && *q != '-') return -2;
        oriv = (*q != '+');
      } else if (i == 2) {
        segw = q;
      } else if (i == 3) {
        if (*q != '+' && *q != '-') return -2;
        oriw = (*q != '+');
      } else if (i == 4) {
        if (*q == '*') {
          ov = ow = 0;
        } else if (*q == ':') {
          ov = INT32_MAX;
          ow = isdigit(*(q + 1)) ? strtol(q + 1, &q, 10) : INT32_MAX;
        } else if (isdigit(*q)) {
          char *r;
          ov = strtol(q, &r, 10);
          if (isupper(*r)) { // CIGAR
            ov = ow = 0;
            do {
              long l;
              l = strtol(q, &q, 10);
              if (*q == 'M' || *q == 'D' || *q == 'N') ov += l;
              if (*q == 'M' || *q == 'I' || *q == 'S') ow += l;
              ++q;
            } while (isdigit(*q));
          } else if (*r == ':') { // overlap lengths
            ow = isdigit(*(r + 1)) ? strtol(r + 1, &r, 10) : INT32_MAX;
          } else break;
        } else break;
        is_ok = 1, rest = c ? p + 1 : 0;
        break;
      }
      ++i, q = p + 1;
      if (c == 0) break;
    }
  }
  if (i == 4 && is_ok == 0) ov = ow = 0, is_ok = 1; // no overlap field
  if (is_ok) {
    uint32_t v, w;
    int l_aux, m_aux = 0;
    uint8_t *aux = 0;
    gfa_arc_t *arc;
    v = gfa_add_seg(g, segv) << 1 | oriv;
    w = gfa_add_seg(g, segw) << 1 | oriw;
    arc = gfa_add_arc1(g, v, w, ov, ow, -1, 0);
    l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
    if (l_aux) {
      gfa_aux_t *a = &g->link_aux[arc->link_id];
      uint8_t *s_L1, *s_L2, *s_SR;
      a->l_aux = l_aux, a->m_aux = m_aux, a->aux = aux;
      s_SR = gfa_aux_get(a->l_aux, a->aux, "SR");
      if (s_SR && s_SR[0] == 'i') {
        arc->rank = *(int32_t *) (s_SR + 1);
        a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_SR);
      }
      s_L1 = gfa_aux_get(a->l_aux, a->aux, "L1");
      if (s_L1) {
        if (ov != INT32_MAX && s_L1[0] == 'i')
          g->seg[v >> 1].len =
                  g->seg[v >> 1].len > ov + *(int32_t *) (s_L1 + 1) ? g->seg[v >> 1].len : ov + *(int32_t *) (s_L1 + 1);
        a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_L1);
      }
      s_L2 = gfa_aux_get(a->l_aux, a->aux, "L2");
      if (s_L2) {
        if (ow != INT32_MAX && s_L2[0] == 'i')
          g->seg[w >> 1].len =
                  g->seg[w >> 1].len > ow + *(int32_t *) (s_L2 + 1) ? g->seg[w >> 1].len : ow + *(int32_t *) (s_L2 + 1);
        a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_L2);
      }
#ifdef GFA_EXT_L
      uint8_t *s_EX = gfa_aux_get(a->l_aux, a->aux, "EX");
      if (s_EX) {
        arc->genes = realloc(arc->genes, strlen((char *)s_EX+1));
        strcpy(arc->genes, (char *) s_EX+1);
        a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_EX);

        /* Alternatively we could keep the string in a->aux 
         * and use a pointer to it in the struct: 
         * 
         * gfa.h: 
         *      gfa_arc_t::exp = unint8_t* 
         * 
         * and here:
         *      arc->exp = s_EX+1;
         * 
         * This way we keep the information inside of `gfa_aux_t`
         * and not in `gfa_arc_t`, so less memory in the arc,
         * but we pay in accessing information.
         */

      } else {
        arc->genes = '\0';
      }
#endif
      if (a->l_aux == 0) {
        free(a->aux);
        a->aux = 0, a->m_aux = 0;
      }
    }
  } else return -1;
  return 0;
}

int gfa_parse_W(gfa_t *g, char *s);
