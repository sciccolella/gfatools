#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include "kstring.h"
#include "gfa-priv.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 65536)

#pragma weak gfa_parse_S
#pragma weak gfa_parse_L
#pragma weak gfa_parse_W

/***********
 * Tag I/O *
 ***********/

int gfa_aux_parse(char *s, uint8_t **data, int *max)
{
	char *q, *p;
	kstring_t str;
	if (s == 0) return 0;
	str.l = 0, str.m = *max, str.s = (char*)*data;
	if (*s == '\t') ++s;
	for (p = q = s;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (p - q >= 5 && q[2] == ':' && q[4] == ':' && (q[3] == 'A' || q[3] == 'i' || q[3] == 'f' || q[3] == 'Z' || q[3] == 'B')) {
				int type = q[3];
				kputsn_(q, 2, &str);
				q += 5;
				if (type == 'A') {
					kputc_('A', &str);
					kputc_(*q, &str);
				} else if (type == 'i') {
					int32_t x;
					x = strtol(q, &q, 10);
					kputc_(type, &str); kputsn_((char*)&x, 4, &str);
				} else if (type == 'f') {
					float x;
					x = strtod(q, &q);
					kputc_('f', &str); kputsn_(&x, 4, &str);
				} else if (type == 'Z') {
					kputc_('Z', &str); kputsn_(q, p - q + 1, &str); // note that this include the trailing NULL
				} else if (type == 'B') {
					type = *q++; // q points to the first ',' following the typing byte
					if (p - q >= 2 && (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I' || type != 'f')) {
						int32_t n;
						char *r;
						for (r = q, n = 0; *r; ++r)
							if (*r == ',') ++n;
						kputc_('B', &str); kputc_(type, &str); kputsn_(&n, 4, &str);
						// TODO: to evaluate which is faster: a) aligned array and then memmove(); b) unaligned array; c) kputsn_()
						if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						else if (type == 'S') while (q + 1 < p) { uint16_t x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						else if (type == 'i') while (q + 1 < p) { int32_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'I') while (q + 1 < p) { uint32_t x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'f') while (q + 1 < p) { float    x = strtod(q + 1, &q);    kputsn_(&x, 4, &str); }
					}
				} // should not be here, as we have tested all types
			}
			q = p + 1;
			if (c == 0) break;
		}
	}
	if (str.l > 0 && str.l == str.m) ks_resize(&str, str.l + 1);
	if (str.s) str.s[str.l] = 0;
	*max = str.m, *data = (uint8_t*)str.s;
	return str.l;
}

int gfa_aux_format(int l_aux, const uint8_t *aux, char **t, int *max)
{
	kstring_t str;
	const uint8_t *s = aux;
	str.l = 0, str.s = *t, str.m = *max;
	while (s < aux + l_aux) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s++;
		kputc('\t', &str); kputsn((char*)key, 2, &str); kputc(':', &str);
		if (type == 'A') { kputsn("A:", 2, &str); kputc(*s, &str); ++s; }
		else if (type == 'i') { kputsn("i:", 2, &str); kputw(*(int32_t*)s, &str); s += 4; }
		else if (type == 'f') { ksprintf(&str, "f:%g", *(float*)s); s += 4; }
		else if (type == 'Z') { kputc(type, &str); kputc(':', &str); while (*s) kputc(*s++, &str); ++s; }
		else if (type == 'B') {
			uint8_t sub_type = *(s++);
			int32_t i, n;
			memcpy(&n, s, 4);
			s += 4; // no point to the start of the array
			kputsn("B:", 2, &str); kputc(sub_type, &str); // write the typing
			for (i = 0; i < n; ++i) { // FIXME: for better performance, put the loop after "if"
				kputc(',', &str);
				if ('c' == sub_type)      { kputw(*(int8_t*)s, &str); ++s; }
				else if ('C' == sub_type) { kputw(*(uint8_t*)s, &str); ++s; }
				else if ('s' == sub_type) { kputw(*(int16_t*)s, &str); s += 2; }
				else if ('S' == sub_type) { kputw(*(uint16_t*)s, &str); s += 2; }
				else if ('i' == sub_type) { kputw(*(int32_t*)s, &str); s += 4; }
				else if ('I' == sub_type) { kputuw(*(uint32_t*)s, &str); s += 4; }
				else if ('f' == sub_type) { ksprintf(&str, "%g", *(float*)s); s += 4; }
			}
		}
	}
	*t = str.s, *max = str.m;
	return str.l;
}

/****************
 * Line parsers *
 ****************/

int gfa_parse_S(gfa_t *g, char *s)
{
	int i, is_ok = 0;
	char *p, *q, *seg = 0, *seq = 0, *rest = 0;
	uint32_t sid, len = 0;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) seg = q;
			else if (i == 1) {
				seq = q[0] == '*'? 0 : gfa_strdup(q);
				is_ok = 1, rest = c? p + 1 : 0;
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
		s_LN = l_aux? gfa_aux_get(l_aux, aux, "LN") : 0;
		if (s_LN && s_LN[0] == 'i') {
			LN = *(int32_t*)(s_LN + 1);
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
				s->snid = gfa_sseq_add(g, (char*)(s_SN + 1)), s->soff = 0;
				l_aux = gfa_aux_del(l_aux, aux, s_SN);
				s_SO = gfa_aux_get(l_aux, aux, "SO");
				if (s_SO && *s_SO == 'i') {
					s->soff = *(int32_t*)(s_SO + 1);
					l_aux = gfa_aux_del(l_aux, aux, s_SO);
				}
			}
			s_SR = gfa_aux_get(l_aux, aux, "SR");
			if (s_SR && *s_SR == 'i') {
				s->rank = *(int32_t*)(s_SR + 1);
				if (s->rank > g->max_rank) g->max_rank = s->rank;
				l_aux = gfa_aux_del(l_aux, aux, s_SR);
			}
			gfa_sseq_update(g, s);
		}
		if (l_aux > 0)
			s->aux.m_aux = m_aux, s->aux.l_aux = l_aux, s->aux.aux = aux;
		else if (aux)
			free(aux);
	} else return -1;
	return 0;
}

int gfa_parse_L(gfa_t *g, char *s)
{
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
					ow = isdigit(*(q+1))? strtol(q+1, &q, 10) : INT32_MAX;
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
						ow = isdigit(*(r+1))? strtol(r+1, &r, 10) : INT32_MAX;
					} else break;
				} else break;
				is_ok = 1, rest = c? p + 1 : 0;
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
				arc->rank = *(int32_t*)(s_SR+1);
				a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_SR);
			}
			s_L1 = gfa_aux_get(a->l_aux, a->aux, "L1");
			if (s_L1) {
				if (ov != INT32_MAX && s_L1[0] == 'i')
					g->seg[v>>1].len = g->seg[v>>1].len > ov + *(int32_t*)(s_L1+1)? g->seg[v>>1].len : ov + *(int32_t*)(s_L1+1);
				a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_L1);
			}
			s_L2 = gfa_aux_get(a->l_aux, a->aux, "L2");
			if (s_L2) {
				if (ow != INT32_MAX && s_L2[0] == 'i')
					g->seg[w>>1].len = g->seg[w>>1].len > ow + *(int32_t*)(s_L2+1)? g->seg[w>>1].len : ow + *(int32_t*)(s_L2+1);
				a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_L2);
			}
			if (a->l_aux == 0) {
				free(a->aux);
				a->aux = 0, a->m_aux = 0;
			}
		}
	} else return -1;
	return 0;
}

int gfa_parse_W(gfa_t *g, char *s)
{
	char *p, *q, *ctg = 0, *sample = 0;
	int32_t i, is_ok = 0;
	char *rest;
	gfa_walk_t t;
	GFA_BZERO(&t, 1);
	for (p = q = s + 2, i = 0;; ++p) {
		t.sample = 0;
		if (*p == 0 || *p == '\t') {
			int32_t c = *p;
			*p = 0;
			if (i == 0) {
				sample = q;
			} else if (i == 1) {
				t.hap = atoi(q);
			} else if (i == 2) {
				ctg = q;
			} else if (i == 3) {
				t.st = atol(q);
			} else if (i == 4) {
				t.en = atol(q);
			} else if (i == 5) {
				char *pp, *qq;
				for (pp = q, t.n_v = 0; pp < p; ++pp)
					if (*pp == '>' || *pp == '<')
						t.n_v++;
				GFA_MALLOC(t.v, t.n_v);
				for (qq = q, pp = q + 1, t.n_v = 0; pp <= p; ++pp) {
					if (pp == p || *pp == '>' || *pp == '<') {
						int32_t a = *pp, seg;
						*pp = 0;
						seg = gfa_name2id(g, qq + 1);
						if (seg >= 0) {
							t.v[t.n_v++] = (uint32_t)seg<<1 | (*qq == '<');
						} else {
							if (gfa_verbose >= 2)
								fprintf(stderr, "WARNING: failed to find segment '%s'\n", qq + 1);
						}
						*pp = a, qq = pp;
					}
				}
				is_ok = 1, rest = c? p + 1 : 0;
				break;
			}
			q = p + 1, ++i;
			if (c == 0) break;
		}
	}
	if (is_ok) {
		int l_aux, m_aux = 0;
		uint8_t *aux = 0;
		l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
		t.sample = gfa_sample_add(g, sample);
		t.snid = gfa_sseq_add(g, ctg);
		if (l_aux > 0)
			t.aux.m_aux = m_aux, t.aux.l_aux = l_aux, t.aux.aux = aux;
		else if (aux)
			free(aux);
		GFA_GROW(gfa_walk_t, g->walk, g->n_walk, g->m_walk);
		g->walk[g->n_walk++] = t;
	} else return -1;
	return 0;
}

static gfa_seg_t *gfa_parse_fa_hdr(gfa_t *g, char *s)
{
	int32_t i;
	char buf[16];
	gfa_seg_t *seg;
	for (i = 0; s[i]; ++i)
		if (isspace(s[i])) break;
	s[i] = 0;
	sprintf(buf, "s%d", g->n_seg + 1);
	i = gfa_add_seg(g, buf);
	seg = &g->seg[i];
	seg->snid = gfa_sseq_add(g, s + 1);
	seg->soff = seg->rank = 0;
	return seg;
}

static void gfa_update_fa_seq(gfa_t *g, gfa_seg_t *seg, int32_t l_seq, const char *seq)
{
	if (seg == 0) return;
	seg->seq = gfa_strdup(seq);
	seg->len = l_seq;
	gfa_sseq_update(g, seg);
}

/****************
 * User-end I/O *
 ****************/

gfa_t *gfa_read(const char *fn)
{
	gzFile fp;
	gfa_t *g;
	kstring_t s = {0,0,0}, fa_seq = {0,0,0};
	kstream_t *ks;
	int dret, is_fa = 0;
	gfa_seg_t *fa_seg = 0;
	uint64_t lineno = 0;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	g = gfa_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
		int ret = 0;
		++lineno;
		if (s.l > 0 && s.s[0] == '>') { // FASTA header
			is_fa = 1;
			if (fa_seg) gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s);
			fa_seg = gfa_parse_fa_hdr(g, s.s);
			fa_seq.l = 0;
		} else if (is_fa) { // FASTA mode
			if (s.l >= 3 && s.s[1] == '\t') { // likely a GFA line
				gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s); // finalize fa_seg
				fa_seg = 0;
				is_fa = 0;
			} else kputsn(s.s, s.l, &fa_seq); // likely a FASTA sequence line
		}
		if (is_fa) continue;
		if (s.l < 3 || s.s[1] != '\t') continue; // empty line
		if (s.s[0] == 'S') ret = gfa_parse_S(g, s.s);
		else if (s.s[0] == 'L') ret = gfa_parse_L(g, s.s);
		else if (s.s[0] == 'W') ret = gfa_parse_W(g, s.s);
		if (ret < 0 && gfa_verbose >= 1)
			fprintf(stderr, "[E] invalid %c-line at line %ld (error code %d)\n", s.s[0], (long)lineno, ret);
	}
	if (is_fa && fa_seg) gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s);
	free(fa_seq.s);
	free(s.s);
	gfa_finalize(g);
	ks_destroy(ks);
	gzclose(fp);
	return g;
}

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		GFA_REALLOC(s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

void gfa_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[32]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 'l' && *(p+1) == 'd') {
				int c, i, l = 0;
				unsigned long x;
				c = va_arg(ap, long);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				++p;
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else {
				fprintf(stderr, "ERROR: unrecognized type '%%%c'\n", *p);
				abort();
			}
			q = p + 1;
		}
	}
	if (p > q) str_copy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}

static void gfa_write_S(kstring_t *out, const gfa_t *g, const gfa_seg_t *s, int flag)
{
	if (s->del) return;
	gfa_sprintf_lite(out, "S\t%s\t", s->name);
	if (s->seq && !(flag & GFA_O_NO_SEQ))
		gfa_sprintf_lite(out, "%s", s->seq);
	else gfa_sprintf_lite(out, "*");
	gfa_sprintf_lite(out, "\tLN:i:%d", s->len);
	if (s->snid >= 0 && s->soff >= 0)
		gfa_sprintf_lite(out, "\tSN:Z:%s\tSO:i:%d", g->sseq[s->snid].name, s->soff);
	if (s->rank >= 0)
		gfa_sprintf_lite(out, "\tSR:i:%d", s->rank);
	if (s->utg && s->utg->n) gfa_sprintf_lite(out, "\tRC:i:%d\tlc:i:%d", s->utg->n, s->utg->len_comp);
	if (s->aux.l_aux > 0) {
		char *t = 0;
		int max = 0;
		gfa_aux_format(s->aux.l_aux, s->aux.aux, &t, &max);
		gfa_sprintf_lite(out, "%s", t);
		free(t);
	}
	gfa_sprintf_lite(out, "\n");
	if (s->utg && s->utg->n) {
		uint32_t j, l;
		for (j = l = 0; j < s->utg->n; ++j) {
			const gfa_utg_t *u = s->utg;
			gfa_sprintf_lite(out, "A\t%s\t%d\t%c\t%s\t%d\t%d\n", s->name, l, "+-"[u->a[j]>>32&1], u->name[j], (int32_t)(u->r[j]>>32), (int32_t)u->r[j]);
			l += (uint32_t)u->a[j];
		}
	}
}

static void gfa_write_L(kstring_t *out, const gfa_t *g, const gfa_arc_t *a, int flag)
{
	const gfa_aux_t *aux = a->link_id < g->n_arc? &g->link_aux[a->link_id] : 0;
	if (a->del || a->comp) return;
	gfa_sprintf_lite(out, "L\t%s\t%c\t%s\t%c", g->seg[a->v_lv>>33].name, "+-"[a->v_lv>>32&1], g->seg[a->w>>1].name, "+-"[a->w&1]);
	if (!(flag & GFA_O_OV_EXT)) {
		gfa_sprintf_lite(out, "\t%dM", a->ov < a->ow? a->ov : a->ow);
	} else {
		if (a->ov == a->ow) gfa_sprintf_lite(out, "\t%dM", a->ov);
		else gfa_sprintf_lite(out, "\t%d:%d", a->ov, a->ow);
	}
	if (a->rank >= 0) gfa_sprintf_lite(out, "\tSR:i:%d", a->rank);
	gfa_sprintf_lite(out, "\tL1:i:%d", gfa_arc_len(*a));
	gfa_sprintf_lite(out, "\tL2:i:%d", gfa_arc_lw(g, *a));
	if (aux && aux->l_aux) {
		char *t = 0;
		int max = 0;
		gfa_aux_format(aux->l_aux, aux->aux, &t, &max);
		if (t) gfa_sprintf_lite(out, "%s", t);
		free(t);
	}
	gfa_sprintf_lite(out, "\n");
}

static void gfa_write_W(kstring_t *out, const gfa_t *g, const gfa_walk_t *w, int flag)
{
	int32_t j;
	gfa_sprintf_lite(out, "W\t%s\t%d\t%s\t", w->sample, w->hap, g->sseq[w->snid].name);
	if (w->st >= 0 && w->en >= 0) gfa_sprintf_lite(out, "%ld\t%ld\t", (long)w->st, (long)w->en);
	else gfa_sprintf_lite(out, "*\t*\t");
	for (j = 0; j < w->n_v; ++j)
		gfa_sprintf_lite(out, "%c%s", "><"[w->v[j]&1], g->seg[w->v[j]>>1].name);
	if (w->aux.l_aux > 0) {
		char *t = 0;
		int max = 0;
		gfa_aux_format(w->aux.l_aux, w->aux.aux, &t, &max);
		gfa_sprintf_lite(out, "%s", t);
		free(t);
	}
	gfa_sprintf_lite(out, "\n");
}

char *gfa_write(const gfa_t *g, int flag, int *len)
{
	uint32_t i;
	uint64_t k;
	kstring_t out = {0,0,0};
	for (i = 0; i < g->n_seg; ++i)
		gfa_write_S(&out, g, &g->seg[i], flag);
	for (k = 0; k < g->n_arc; ++k)
		gfa_write_L(&out, g, &g->arc[k], flag);
	for (i = 0; i < g->n_walk; ++i)
		gfa_write_W(&out, g, &g->walk[i], flag);
	*len = out.l;
	return out.s;
}

void gfa_print(const gfa_t *g, FILE *fp, int flag)
{
	uint32_t i;
	uint64_t k;
	kstring_t out = {0,0,0};
	// S-lines
	for (i = 0; i < g->n_seg; ++i) {
		out.l = 0;
		gfa_write_S(&out, g, &g->seg[i], flag);
		fwrite(out.s, 1, out.l, fp);
	}
	// L-lines
	for (k = 0; k < g->n_arc; ++k) {
		out.l = 0;
		gfa_write_L(&out, g, &g->arc[k], flag);
		fwrite(out.s, 1, out.l, fp);
	}
	// W-lines
	for (i = 0; i < g->n_walk; ++i) {
		out.l = 0;
		gfa_write_W(&out, g, &g->walk[i], flag);
		fwrite(out.s, 1, out.l, fp);
	}
	free(out.s);
}

/*
 * Read a list of labels
 */
char **gfa_read_list(const char *o, int *n_)
{
	int n = 0, m = 0;
	char **s = 0;
	*n_ = 0;
	if (*o != '@') {
		const char *q = o, *p;
		for (p = q;; ++p) {
			if (*p == ',' || *p == ' ' || *p == '\t' || *p == 0) {
				if (p - q > 0) {
					GFA_GROW0(char*, s, n, m);
					s[n++] = gfa_strndup(q, p - q);
				}
				if (*p == 0) break;
				q = p + 1;
			}
		}
	} else {
		gzFile fp;
		kstream_t *ks;
		kstring_t str = {0,0,0};
		int dret;

		fp = gzopen(o + 1, "r");
		if (fp == 0) return 0;
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
			char *p;
			for (p = str.s; *p && !isspace(*p); ++p);
			GFA_GROW0(char*, s, n, m);
			s[n++] = gfa_strndup(str.s, p - str.s);
		}
		ks_destroy(ks);
		gzclose(fp);
	}
	if (s) s = (char**)realloc(s, n * sizeof(char*));
	*n_ = n;
	return s;
}
