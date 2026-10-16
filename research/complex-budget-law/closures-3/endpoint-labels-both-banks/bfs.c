// expQW bfs.c: exact minimum loss L of a tiny both-bank one-block CM instance (see mkinst.py), by bucketed
// Dijkstra over (frames, value matrix over GF(P)). Prints min L found <= LMAX, or "none <= LMAX".
// usage: bfs inst.txt LMAX log2cap
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define MAXF 64
#define MAXN 6
int nf, n, P, su, dimf[MAXF]; char sub[MAXF][MAXF]; int st[MAXN], sk[MAXN], tc[MAXN], tv[MAXN];
typedef struct { uint8_t fr[MAXN]; uint8_t m[MAXN][MAXN]; } S;
static uint64_t enc(const S *s) { uint64_t k = 0; for (int i = 0; i < n; i++) k = k * nf + s->fr[i];
  for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) k = k * P + s->m[i][j]; return k; }
static void dec(uint64_t k, S *s) { for (int i = n - 1; i >= 0; i--) for (int j = n - 1; j >= 0; j--) { s->m[i][j] = k % P; k /= P; }
  for (int i = n - 1; i >= 0; i--) { s->fr[i] = k % nf; k /= nf; } }
uint64_t *tab; uint64_t cap, used = 0;           // entries: key+1 in low 58 bits? store key and cost separately
uint8_t *cost;
static int64_t look(uint64_t key, int c) {     // insert if new or better; return 1 if improved
  uint64_t h = (key * 0x9E3779B97F4A7C15ULL) & (cap - 1);
  while (tab[h]) { if (tab[h] == key + 1) { if (c < cost[h]) { cost[h] = c; return 1; } return 0; } h = (h + 1) & (cap - 1); }
  if (++used > cap * 0.9) { fprintf(stderr, "table full\n"); exit(2); }
  tab[h] = key + 1; cost[h] = c; return 1; }
typedef struct { uint64_t *a; size_t len, cap; } Q;
static void push(Q *q, uint64_t k) { if (q->len == q->cap) { q->cap = q->cap ? 2 * q->cap : 1 << 20; q->a = realloc(q->a, q->cap * 8); } q->a[q->len++] = k; }
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); int LMAX = atoi(argv[2]); int lc = atoi(argv[3]);
  if (fscanf(f, "%d %d %d %d", &nf, &n, &P, &su) != 4) return 1;
  char buf[256]; for (int i = 0; i < nf; i++) { if (fscanf(f, "%d %s", &dimf[i], buf) != 2) return 1; for (int j = 0; j < nf; j++) sub[i][j] = buf[j] == '1'; }
  for (int i = 0; i < n; i++) if (fscanf(f, "%d", &st[i]) != 1) return 1;
  for (int i = 0; i < n; i++) if (fscanf(f, "%d", &sk[i]) != 1) return 1;
  for (int i = 0; i < n; i++) if (fscanf(f, "%d %d", &tc[i], &tv[i]) != 2) return 1;
  cap = 1ULL << lc; tab = calloc(cap, 8); cost = malloc(cap);
  S s0; memset(&s0, 0, sizeof s0); for (int i = 0; i < n; i++) { s0.fr[i] = st[i]; s0.m[i][i] = 1; }
  Q *qs = calloc(LMAX + 1, sizeof(Q)); uint64_t k0 = enc(&s0); look(k0, 0); push(&qs[0], k0);
  for (int L = 0; L <= LMAX; L++) {
    Q *q = &qs[L];
    for (size_t qi = 0; qi < q->len; qi++) {
      uint64_t k = q->a[qi]; S s; dec(k, &s);
      { uint64_t h = (k * 0x9E3779B97F4A7C15ULL) & (cap - 1); while (tab[h] != k + 1) h = (h + 1) & (cap - 1); if (cost[h] < L) continue; }
      int good = 1;                                   // target?
      for (int i = 0; i < n && good; i++) { if (!sub[s.fr[i]][sk[i]]) good = 0;
        for (int j = 0; j < n && good; j++) if (s.m[i][j] != (j == tc[i] ? tv[i] : 0)) good = 0; }
      if (good) { printf("min L = %d  (sum u = %d, states %llu)\n", L, su, (unsigned long long)used); return 0; }
      for (int r = 0; r < n; r++) for (int g = 0; g < nf; g++) {           // frame steps
        int a = s.fr[r]; if (g == a) continue; int c;
        if (sub[a][g]) c = 0; else if (sub[g][a]) c = dimf[a] - dimf[g]; else continue;
        if (L + c > LMAX) continue;
        S t = s; t.fr[r] = g; uint64_t kk = enc(&t); if (look(kk, L + c)) push(&qs[L + c], kk); }
      for (int r = 0; r < n; r++) for (int u = 0; u < n; u++) if (u != r && s.fr[r] == s.fr[u])  // transvections
        for (int c = 1; c < P; c++) { S t = s; for (int j = 0; j < n; j++) t.m[r][j] = (t.m[r][j] + c * t.m[u][j]) % P;
          uint64_t kk = enc(&t); if (look(kk, L)) push(&qs[L], kk); }
    }
    fprintf(stderr, "  level %d done: %zu queued, %llu states\n", L, q->len, (unsigned long long)used);
    free(q->a); q->a = 0;
  }
  printf("none <= %d  (sum u = %d, states %llu)\n", LMAX, su, (unsigned long long)used); return 0;
}
