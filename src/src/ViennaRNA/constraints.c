/* constraints handling */

#include <assert.h>
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/energy_const.h" /* defines MINPSCORE */
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void
hc_add_up(vrna_fold_compound_t *vc,
          int i,
          char option);

PRIVATE INLINE  void
hc_cant_pair( unsigned int i,
              char c_option,
              char *hc,
              unsigned int length,
              unsigned int min_loop_size,
              int *index);

PRIVATE INLINE  void
hc_must_pair( unsigned int i,
              char c_option,
              char *hc,
              int *index);

PRIVATE INLINE  void
hc_pairs_upstream(unsigned int i,
                  char c_option,
                  char *hc,
                  unsigned int length,
                  int *index);

PRIVATE INLINE  void
hc_pairs_downstream(unsigned int i,
                    char c_option,
                    char *hc,
                    unsigned int length,
                    int *index);

PRIVATE INLINE  void
hc_allow_pair(unsigned int i,
              unsigned int j,
              char c_option,
              char *hc,
              int *index);

PRIVATE INLINE  void
hc_weak_enforce_pair( unsigned int i,
                      unsigned int j,
                      char c_option,
                      char *hc,
                      unsigned int length,
                      unsigned int min_loop_size,
                      int *index);

PRIVATE INLINE  void
hc_enforce_pair(unsigned int i,
                unsigned int j,
                char c_option,
                char *hc,
                unsigned int length,
                unsigned int min_loop_size,
                int *index);

PRIVATE INLINE  void
hc_intramolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index);

PRIVATE INLINE  void
hc_intermolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index);

PRIVATE INLINE  void
adjust_ptypes(char *ptype,
              vrna_hc_t *hc,
              unsigned int length,
              unsigned int indx_type);

PRIVATE void
apply_DB_constraint(const char *constraint,
                    char *ptype,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int cut,
                    unsigned int options);

PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc);

PRIVATE void
hc_update_up(vrna_fold_compound_t *vc);

PRIVATE void
sc_parse_parameters(const char *string,
                    char c1,
                    char c2,
                    float *v1,
                    float *v2);

PRIVATE void
sc_add_up_mfe(vrna_fold_compound_t *vc,
              const FLT_OR_DBL *constraints,
              unsigned int options);

PRIVATE void
sc_add_up_pf( vrna_fold_compound_t *vc,
              const FLT_OR_DBL *constraints,
              unsigned int options);

PRIVATE void
sc_add_bp_mfe(vrna_fold_compound_t *vc,
              const FLT_OR_DBL **constraints,
              unsigned int options);

PRIVATE void
sc_add_bp_pf( vrna_fold_compound_t *vc,
              const FLT_OR_DBL **constraints,
              unsigned int options);

PRIVATE void
sc_add_stack_en_mfe(vrna_fold_compound_t *vc,
                    const FLT_OR_DBL *constraints,
                    unsigned int options);

PRIVATE void
sc_add_stack_en_pf( vrna_fold_compound_t *vc,
                    const FLT_OR_DBL *constraints,
                    unsigned int options);


/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
PUBLIC  void
vrna_message_constraint_options_all(void){

  vrna_message_constraint_options(  VRNA_CONSTRAINT_DB_PIPE
                                  | VRNA_CONSTRAINT_DB_DOT
                                  | VRNA_CONSTRAINT_DB_X
                                  | VRNA_CONSTRAINT_DB_ANG_BRACK
                                  | VRNA_CONSTRAINT_DB_RND_BRACK);
}

PUBLIC  void
vrna_message_constraint_options(unsigned int option){

  if(!(option & VRNA_CONSTRAINT_NO_HEADER)) printf("Input structure constraints using the following notation:\n");
  if(option & VRNA_CONSTRAINT_DB_PIPE)       printf("| : paired with another base\n");
  if(option & VRNA_CONSTRAINT_DB_DOT)        printf(". : no constraint at all\n");
  if(option & VRNA_CONSTRAINT_DB_X)          printf("x : base must not pair\n");
  if(option & VRNA_CONSTRAINT_DB_ANG_BRACK)  printf("< : base i is paired with a base j<i\n> : base i is paired with a base j>i\n");
  if(option & VRNA_CONSTRAINT_DB_RND_BRACK)  printf("matching brackets ( ): base i pairs base j\n");
}

PUBLIC  void
vrna_constraints_add( vrna_fold_compound_t *vc,
                      const char *constraint,
                      unsigned int options){

  int         i, d;
  vrna_md_t   *md;

  if(vc){
    if(vc->params)
      md = &(vc->params->model_details);
    else if(vc->exp_params)
      md = &(vc->exp_params->model_details);
    else
      vrna_message_error("constraints.c@vrna_add_constraints: fold compound has no params or exp_params");

    if(!vc->hc)
      vrna_hc_init(vc);

    if(options & VRNA_CONSTRAINT_DB){ /* apply hard constraints from dot-bracket notation */
      apply_DB_constraint(constraint,
                          vc->hc->matrix,
                          vc->length,
                          md->min_loop_size,
                          -1,
                          options);
      hc_update_up(vc);
    } else if(options & VRNA_CONSTRAINT_FILE){ /* constraints from file */
      plist *p, *c = vrna_file_constraints_read(constraint, vc->length, options);

      /* now do something with the constraints we've just read */
      if(c){
        FLT_OR_DBL  **sc_bp       = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (vc->length + 1));
        FLT_OR_DBL  *sc_up        = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));
        int     sc_up_present = 0;
        int     sc_bp_present = 0;

        for(i = 0; i <= vc->length; i++)
          sc_bp[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));

        for(p = c; p->i; p++){
          if(p->type & 4096){ /* soft constraint */
            if(p->j == 0){  /* pseudo energy for unpairedness */
              sc_up_present = 1;
              sc_up[p->i] += (FLT_OR_DBL)p->p;
            } else { /* pseudo energy for base pair */
              sc_bp_present = 1;
              sc_bp[p->i][p->j] += (FLT_OR_DBL)p->p;
            }
          } else {  /* hard constraint */
            if(p->j == 0){
              hc_add_up(vc, p->i, (char)p->type);
            } else if(p->i == p->j){ 
              d = 0;
              if(1024 & p->type)
                d = -1;
              else if(2048 & p->type)
                d = 1;
              vrna_hc_add_bp_nonspecific(vc, p->i, d, (char)(p->type));
            } else {
              if(p->type & 8192){
                vrna_hc_add_bp(vc, p->i, p->j, (char)(p->type) | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE);
              } else {
                vrna_hc_add_bp(vc, p->i, p->j, (char)(p->type));
              }
            }
          }
        }

        hc_update_up(vc);

        /* ############################### */
        /* init empty soft constraints     */
        /* ############################### */
        if(sc_up_present || sc_bp_present){
          vrna_sc_init(vc);
          if(sc_bp_present)
            vrna_sc_add_bp(vc, (const FLT_OR_DBL **)sc_bp, options);
          if(sc_up_present)
            vrna_sc_add_up(vc, (const FLT_OR_DBL *)sc_up, options);
        }
        /* clean up */
        for(i = 0; i <= vc->length; i++)
          free(sc_bp[i]);
        free(sc_bp);
        free(sc_up);
      }

      free(c);
    }
  }
}

PUBLIC  void
vrna_hc_init(vrna_fold_compound_t *vc){

  unsigned int  n;
  vrna_hc_t     *hc;

  n           = vc->length;

  /* free previous hard constraints */
  vrna_hc_free(vc->hc);

  /* allocate memory new hard constraints data structure */
  hc          = (vrna_hc_t *)vrna_alloc(sizeof(vrna_hc_t));
  hc->matrix  = (char *)vrna_alloc(sizeof(char)*((n*(n+1))/2+2));
  hc->up_ext  = (int *)vrna_alloc(sizeof(int)*(n+2));
  hc->up_hp   = (int *)vrna_alloc(sizeof(int)*(n+2));
  hc->up_int  = (int *)vrna_alloc(sizeof(int)*(n+2));
  hc->up_ml   = (int *)vrna_alloc(sizeof(int)*(n+2));

  /* set new hard constraints */
  vc->hc = hc;

  /* prefill default values  */
  hc_reset_to_default(vc);

  /* add null pointers for the generalized hard constraint feature */
  hc->f           = NULL;
  hc->data        = NULL;
  hc->free_data   = NULL;

  /* update */
  hc_update_up(vc);
}

PUBLIC void
vrna_hc_add_up( vrna_fold_compound_t *vc,
                int i,
                char option){

  int j;

  if(vc)
    if(vc->hc){
      if((i <= 0) || (i > vc->length)){
        vrna_message_warning("vrna_hc_add_up: position out of range, not doing anything");
        return;
      }

      hc_add_up(vc, i, option);

      hc_update_up(vc);
    }
}

PRIVATE void
hc_add_up(vrna_fold_compound_t *vc,
          int i,
          char option){

  int   j;
  char  type = (char)0;

  if(option & VRNA_CONSTRAINT_CONTEXT_ENFORCE){ /* force nucleotide to appear unpaired within a certain type of loop */
    /* do not allow i to be paired with any other nucleotide */
    if(!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)){
      for(j = 1; j < i; j++)
        vc->hc->matrix[vc->jindx[i] + j] = (char)0;
      for(j = i+1; j <= vc->length; j++)
        vc->hc->matrix[vc->jindx[j] + i] = (char)0;
    }

    type = option & (char)( VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                            | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                            | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                            | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);

    vc->hc->matrix[vc->jindx[i] + i] = type;
  } else {
    type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

    /* do not allow i to be paired with any other nucleotide (in context type) */
    if(!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)){
      for(j = 1; j < i; j++)
        vc->hc->matrix[vc->jindx[i] + j] &= ~type;
      for(j = i+1; j <= vc->length; j++)
        vc->hc->matrix[vc->jindx[j] + i] &= ~type;
    }

    vc->hc->matrix[vc->jindx[i] + i] = (char)(  VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                                              | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                                              | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                                              | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
  }
}

PUBLIC void
vrna_hc_add_bp_nonspecific( vrna_fold_compound_t *vc,
                            int i,
                            int d,
                            char option){
  int   p;
  char  type, t1, t2;

  if(vc)
    if(vc->hc){
      if((i <= 0) || (i > vc->length)){
        vrna_message_warning("vrna_hc_add_bp_nonspecific: position out of range, not doing anything");
        return;
      }

      /* force position i to pair with some other nucleotide */
      type  = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
      /* force direction */
      t1    = (d <= 0) ? type : (char)0;
      t2    = (d >= 0) ? type : (char)0;
      for(p = 1; p < i; p++)
        vc->hc->matrix[vc->jindx[i] + p] &= t1;
      for(p = i+1; p <= vc->length; p++)
        vc->hc->matrix[vc->jindx[p] + i] &= t2;

      vc->hc->matrix[vc->jindx[i] + i] = (char)0;

      hc_update_up(vc);
    }

}

PUBLIC void
vrna_hc_add_bp( vrna_fold_compound_t *vc,
                int i,
                int j,
                char option){

  int   k, l;
  char  type;

  if(vc)
    if(vc->hc){
      if((i <= 0) || (j <= i) || (j > vc->length)){
        vrna_message_warning("vrna_hc_add_bp: position out of range, not doing anything");
        return;
      }

      /* reset ptype in case (i,j) is a non-canonical pair */
      if(option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS){
        if(vc->hc->matrix[vc->jindx[j] + i])
          if(vc->ptype[vc->jindx[j] + i] == 0)
            vc->ptype[vc->jindx[j] + i] = 7;
      }

      vc->hc->matrix[vc->jindx[j] + i] = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

      if(!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)){
        /*
          remove all conflicting base pairs, i.e. do not allow i,j to pair
          with any other nucleotide k
        */
        for(k = 1; k < i; k++){
          vc->hc->matrix[vc->jindx[i] + k] = (char)0;
          vc->hc->matrix[vc->jindx[j] + k] = (char)0;
          for(l = i+1; l < j; l++)
            vc->hc->matrix[vc->jindx[l] + k] = (char)0;
        }
        for(k = i+1; k < j; k++){
          vc->hc->matrix[vc->jindx[k] + i] = (char)0;
          vc->hc->matrix[vc->jindx[j] + k] = (char)0;
          for(l = j + 1; l <= vc->length; l++)
            vc->hc->matrix[vc->jindx[l] + k] = (char)0;
        }
        for(k = j+1; k <= vc->length; k++){
          vc->hc->matrix[vc->jindx[k] + i] = (char)0;
          vc->hc->matrix[vc->jindx[k] + j] = (char)0;
        }
      }

      if(option & VRNA_CONSTRAINT_CONTEXT_ENFORCE){

        /* do not allow i,j to be unpaired */
        vc->hc->matrix[vc->jindx[i] + i] = (char)0;
        vc->hc->matrix[vc->jindx[j] + j] = (char)0;

        hc_update_up(vc);
      }
    }
}

PUBLIC void
vrna_hc_free(vrna_hc_t *hc){

  if(hc){
    free(hc->matrix);
    free(hc->up_ext);
    free(hc->up_hp);
    free(hc->up_int);
    free(hc->up_ml);

    if(hc->free_data)
      hc->free_data(hc->data);

    free(hc);
  }
}

#ifdef WITH_GEN_HC

PUBLIC void
vrna_hc_add_f(vrna_fold_compound_t *vc,
              vrna_callback_hc_evaluate *f){

  if(vc && f){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->hc)
        vrna_hc_init(vc);

      vc->hc->f = f;
    }
  }
}

PUBLIC void
vrna_hc_add_data( vrna_fold_compound_t *vc,
                  void *data,
                  vrna_callback_free_auxdata *f){

  if(vc && data){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->hc)
        vrna_hc_init(vc);

      vc->hc->data        = data;
      vc->hc->free_data   = f;
    }
  }
}

#endif

PRIVATE void
apply_DB_constraint(const char *constraint,
                    char *hc,
                    unsigned int length,
                    unsigned int min_loop_size,
                    int cut,
                    unsigned int options){

  int n,i,j;
  int hx, *stack;
  int *index;
  char c_option;

  if(constraint == NULL) return;

  n         = (int)strlen(constraint);
  stack     = (int *) vrna_alloc(sizeof(int)*(n+1));
  index     = vrna_idx_col_wise(length);
  c_option  =   VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
              | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
              | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
              | VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC
              | VRNA_CONSTRAINT_CONTEXT_MB_LOOP
              | VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC;

  for(hx=0, j=1; j<=n; j++) {
    switch (constraint[j-1]) {
       /* can't pair */
       case 'x':  if(options & VRNA_CONSTRAINT_DB_X){
                    hc_cant_pair(j, c_option, hc, length, min_loop_size, index);
                  }
                  break;

      /* must pair, i.e. may not be unpaired */
      case '|':   if(options & VRNA_CONSTRAINT_DB_PIPE){
                    if(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
                      hc_must_pair(j, c_option, hc, index);
                  }
                  break;

      /* weak enforced pair 'open' */
      case '(':   if(options & VRNA_CONSTRAINT_DB_RND_BRACK){
                    stack[hx++]=j;
                  }
                  break;

      /* weak enforced pair 'close' */
      case ')':   if(options & VRNA_CONSTRAINT_DB_RND_BRACK){
                    if (hx<=0) {
                      fprintf(stderr, "%s\n", constraint);
                      vrna_message_error("unbalanced brackets in constraints");
                    }
                    i = stack[--hx];
                    if(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
                      hc_enforce_pair(i, j, c_option, hc, length, min_loop_size, index);
                    else
                      hc_weak_enforce_pair(i, j, c_option, hc, length, min_loop_size, index);
                  }
                  break;

      /* pairs upstream */
      case '<':   if(options & VRNA_CONSTRAINT_DB_ANG_BRACK){
                    hc_pairs_downstream(j, c_option, hc, length, index);
                    if(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
                      hc_must_pair(j, c_option, hc, index);
                  }
                  break;

      /* pairs downstream */
      case '>':   if(options & VRNA_CONSTRAINT_DB_ANG_BRACK){
                    hc_pairs_upstream(j, c_option, hc, length, index);
                    if(options & VRNA_CONSTRAINT_DB_ENFORCE_BP)
                      hc_must_pair(j, c_option, hc, index);
                  }
                  break;

      /* only intramolecular basepairing */
      case 'l':   if(options & VRNA_CONSTRAINT_DB_INTRAMOL){
                    hc_intramolecular_only(j, c_option, hc, length, min_loop_size, cut, index);
                  }
                  break;

      /* only intermolecular bp */
      case 'e':   if(options & VRNA_CONSTRAINT_DB_INTERMOL){
                    hc_intermolecular_only(j, c_option, hc, length, min_loop_size, cut, index);
                  }
                  break;

      case '.':   break;

      default:    vrna_message_warning("unrecognized character in pseudo dot-bracket notation constraint string\n");
                  break;
    }
  }

  if (hx!=0) {
    fprintf(stderr, "%s\n", constraint);
    vrna_message_error("unbalanced brackets in constraint string");
  }
  /* clean up */
  free(index);
  free(stack);
}

PRIVATE INLINE  void
hc_intramolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index){

  unsigned int l;

  if(cut > 1){
    if(i < cut)
      for(l = MAX2(i+min_loop_size, cut); l <= length; l++)
        hc[index[l] + i] &= ~c_option;
    else
      for(l = 1; l < MIN2(cut, i-min_loop_size); l++)
        hc[index[i] + l] &= ~c_option;
  }
}

PRIVATE INLINE  void
hc_intermolecular_only( unsigned int i,
                        char c_option,
                        char *hc,
                        unsigned int length,
                        unsigned int min_loop_size,
                        int cut,
                        int *index){

  unsigned int l;

  if(cut > 1){
    if(i < cut){
      for(l = 1; l < i; l++)
        hc[index[i] + l] &= ~c_option;
      for(l = i + 1; l < cut; l++)
        hc[index[l] + i] &= ~c_option;
    } else {
      for(l = cut; l < i; l++)
        hc[index[i] + l] &= ~c_option;
      for(l = i + 1; l <= length; l++)
        hc[index[l] + i] &= ~c_option;
    }
  }
}

PRIVATE INLINE  void
hc_cant_pair( unsigned int i,
              char c_option,
              char *hc,
              unsigned int length,
              unsigned int min_loop_size,
              int *index){

  hc_pairs_upstream(i, c_option, hc, length, index);
  hc_pairs_downstream(i, c_option, hc, length, index);
}

PRIVATE INLINE  void
hc_must_pair( unsigned int i,
              char c_option,
              char *hc,
              int *index){

  hc[index[i]+i] &= ~c_option;
}

PRIVATE INLINE  void
hc_pairs_upstream(unsigned int i,
                  char c_option,
                  char *hc,
                  unsigned int length,
                  int *index){

  unsigned int l;

  /* prohibit downstream pairs */
  for(l = length; l > i; l--)
    hc[index[l] + i] = (char)0;
  /* allow upstream pairs of given type */
  for(l = i - 1; l >= 1; l--)
    hc[index[i] + l] &= c_option;
}

PRIVATE INLINE  void
hc_pairs_downstream(unsigned int i,
                    char c_option,
                    char *hc,
                    unsigned int length,
                    int *index){

  unsigned int l;
  /* allow downstream pairs of given type */
  for(l = length; l > i; l--)
    hc[index[l] + i] &= c_option;
  /* forbid upstream pairs */
  for(l = i - 1; l >= 1; l--)
    hc[index[i] + l] = (char)0;
}

PRIVATE INLINE  void
hc_allow_pair(unsigned int i,
              unsigned int j,
              char c_option,
              char *hc,
              int *index){

  hc[index[j] + i] |= c_option;
}

PRIVATE INLINE  void
hc_weak_enforce_pair( unsigned int i,
                      unsigned int j,
                      char c_option,
                      char *hc,
                      unsigned int length,
                      unsigned int min_loop_size,
                      int *index){

  unsigned int k, l;

  /* don't allow pairs (k,i) 1 <= k < i */
  /* don't allow pairs (i,k) i < k <= n */ 
  hc_pairs_upstream(i, (char)0, hc, length, index);
  /* don't allow pairs (k,j) 1 <= k < j */
  /* don't allow pairs (j,k) j < k <= n */ 
  hc_pairs_upstream(j, (char)0, hc, length, index);

  /* don't allow pairs i < k < j < l */
  for(k = i+1; k < j; k++)
    for(l = j+1; l <= length; l++){
      hc[index[l] + k] = 0;
    }
  /* don't allow pairs k<i<l<j */
  for(k = 1; k < i; k++)
    for(l = i+1; l < j; l++){
      hc[index[l] + k] = 0;
    }
  /* allow base pair (i,j) */
  hc[index[j] + i] |= c_option;
}

PRIVATE INLINE  void
hc_enforce_pair(unsigned int i,
                unsigned int j,
                char c_option,
                char *hc,
                unsigned int length,
                unsigned int min_loop_size,
                int *index){

  hc_weak_enforce_pair( i,
                        j,
                        c_option,
                        hc,
                        length,
                        min_loop_size,
                        index);

  /* forbid i and j to be unpaired */
  hc[index[i] + i] = 0;
  hc[index[j] + j] = 0;
}

PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc){

  unsigned int      i, j, ij, min_loop_size, n;
  int               max_span, *idx;
  vrna_md_t         *md;
  vrna_hc_t         *hc;
  short             *S;

  n   = vc->length;
  hc  = vc->hc;
  idx = vc->jindx;
  S   = vc->sequence_encoding;

  if(vc->params)
    md  = &(vc->params->model_details);
  else if(vc->exp_params)
    md  = &(vc->exp_params->model_details);
  else
    vrna_message_error("missing model_details in fold_compound");

  min_loop_size = md->min_loop_size;
  max_span      = md->max_bp_span;

  if((max_span < 5) || (max_span > n))
    max_span  = n;

  /* ######################### */
  /* fill with default values  */
  /* ######################### */

  /* 1. unpaired nucleotides are allowed in all contexts */
  for(i = 1; i <= n; i++)
    hc->matrix[idx[i] + i]  =   VRNA_CONSTRAINT_CONTEXT_EXT_LOOP
                              | VRNA_CONSTRAINT_CONTEXT_HP_LOOP
                              | VRNA_CONSTRAINT_CONTEXT_INT_LOOP
                              | VRNA_CONSTRAINT_CONTEXT_MB_LOOP;

  if(vc->type == VRNA_VC_TYPE_ALIGNMENT){
    /* 2. all base pairs with pscore above threshold are allowed in all contexts */
    for(j = n; j > min_loop_size + 1; j--){
      ij = idx[j]+1;
      for(i=1; i < j - min_loop_size; i++, ij++)
        if((j-i+1) > max_span){
          hc->matrix[ij] = (char)0;
        } else {
          hc->matrix[ij] = (vc->pscore[idx[j]+i] >= md->cv_fact*MINPSCORE) ? VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS : (char)0;
        }
    }
    /* correct for no lonely pairs (assuming that ptypes already incorporate noLP status) */
    /* this should be included in the pscore which is checked above */
  } else {
    /* 2. all canonical base pairs are allowed in all contexts */
    for(j = n; j > min_loop_size + 1; j--){
      ij = idx[j]+1;
      for(i=1; i < j - min_loop_size; i++, ij++){
        char opt = (char)0;
        if((j-i+1) <= max_span){
          int t = md->pair[S[i]][S[j]];
          switch(t){
            case 0:   break;
            case 3:   /* fallthrough */
            case 4:   if(md->noGU){
                        break;
                      } else if(md->noGUclosure){
                        opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS & ~(VRNA_CONSTRAINT_CONTEXT_HP_LOOP | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
                        break;
                      } /* else fallthrough */
            default:  opt = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
                      break;
          }
        }
        hc->matrix[ij] = opt;
      }
    }

    /* correct for no lonely pairs (assuming that ptypes already incorporate noLP status) */
    /* this should be fixed such that ij loses its hard constraint type if it does not
       allow for enclosing an interior loop, etc.
    */
    if(md->noLP)
      for(i = 1; i < n; i++)
        for(j = i + min_loop_size + 1; j <= n; j++){
          if(hc->matrix[idx[j] +i]){
            if(!vc->ptype[idx[j] + i]){
              hc->matrix[idx[j] + i] = (char)0;
            }
          }
        }
  }

  /* should we reset the generalized hard constraint feature here? */
  if(hc->f || hc->data){
    if(hc->free_data)
      hc->free_data(hc->data);

    hc->f           = NULL;
    hc->data        = NULL;
    hc->free_data   = NULL;
  }

}

PRIVATE void
hc_update_up(vrna_fold_compound_t *vc){

  unsigned int      i, n;
  int               *idx;
  vrna_hc_t         *hc;

  n   = vc->length;
  idx = vc->jindx;
  hc  = vc->hc;

  for(hc->up_ext[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
    hc->up_ext[i] = (hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ? 1 + hc->up_ext[i+1] : 0;

  for(hc->up_hp[n+1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
    hc->up_hp[i] = (hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ? 1 + hc->up_hp[i+1] : 0;

  for(hc->up_int[n+1] = 0, i = n; i > 0; i--) /* unpaired stretch in interior loop */
    hc->up_int[i] = (hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 + hc->up_int[i+1] : 0;

  for(hc->up_ml[n+1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
    hc->up_ml[i] = (hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 + hc->up_ml[i+1] : 0;

  /*
   *  loop arround once more until we find a nucleotide that mustn't
   *  be unpaired (needed for circular folding)
   */

  if(hc->matrix[idx[1]+1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
    hc->up_ext[n+1] = hc->up_ext[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP){
        hc->up_ext[i] = MIN2(n, 1 + hc->up_ext[i+1]);
      } else
        break;
    }
  }

  if(hc->matrix[idx[1]+1] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){
    hc->up_hp[n+1] = hc->up_hp[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP){
        hc->up_hp[i] = MIN2(n, 1 + hc->up_hp[i+1]);
      } else
        break;
    }
  }

  if(hc->matrix[idx[1]+1] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
    hc->up_int[n+1] = hc->up_int[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP){
        hc->up_int[i] = MIN2(n, 1 + hc->up_int[i+1]);
      } else
        break;
    }
  }

  if(hc->matrix[idx[1]+1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
    hc->up_ml[n+1] = hc->up_ml[1];
    for(i = n; i > 0; i--){
      if(hc->matrix[idx[i]+i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP){
        hc->up_ml[i] = MIN2(n, 1 + hc->up_ml[i+1]);
      } else
        break;
    }
  }

}


/* And now, for something completely different...
 *
 * ... the soft constraints section follows below
 */

PUBLIC int
vrna_sc_SHAPE_to_pr(const char *shape_conversion,
                    double *values,
                    int length,
                    double default_value){

  int *indices;
  int i, j;
  int index;
  int ret = 1;

  if(!shape_conversion || !(*shape_conversion) || length <= 0)
    return 0;

  if(*shape_conversion == 'S')
    return 1;

  indices = vrna_alloc(sizeof(int) * (length + 1));
  for (i = 1, j = 0; i <= length; ++i){
    if(values[i] < 0)
      values[i] = default_value;
    else
      indices[j++] = i;
  }

  if(*shape_conversion == 'M'){
    double max;
    double map_info[4][2] = {{0.25, 0.35},
                           {0.30, 0.55},
                           {0.70, 0.85},
                           {0, 1}};

    max = values[1];
    for(i = 2; i <= length; ++i)
      max = MAX2(max, values[i]);
    map_info[3][0] = max;

    for(i = 0; indices[i]; ++i){
      double lower_source = 0;
      double lower_target = 0;

      index = indices[i];

      if(values[index] == 0)
        continue;

      for(j = 0; j < 4; ++j){
        if(values[index] > lower_source && values[index] <= map_info[j][0]){
          double diff_source = map_info[j][0] - lower_source;
          double diff_target = map_info[j][1] - lower_target;
          values[index] = (values[index] - lower_source) / diff_source * diff_target + lower_target;
          break;
        }

        lower_source = map_info[j][0];
        lower_target = map_info[j][1];
      }
    }
  }
  else if (*shape_conversion == 'C'){
    float cutoff = 0.25;
    int i;

    sscanf(shape_conversion + 1, "%f", &cutoff);

    for(i = 0; indices[i]; ++i){
      index = indices[i];
      values[index] = values[index] < cutoff ? 0 : 1;
    }
  }
  else if (*shape_conversion == 'L' || *shape_conversion == 'O'){
    int i;
    float slope = (*shape_conversion == 'L') ? 0.68 : 1.6;
    float intercept = (*shape_conversion == 'L') ? 0.2 : -2.29;

    sc_parse_parameters(shape_conversion + 1, 's', 'i', &slope, &intercept);

    for(i = 0; indices[i]; ++i){
      double v;
      index = indices[i];

      v = (*shape_conversion == 'L') ? values[index] : log(values[index]);
      values[index] = MAX2(MIN2((v - intercept) / slope, 1),0);
    }
  }
  else
    ret = 0;

  free(indices);

  return ret;
}

PUBLIC void
vrna_sc_init(vrna_fold_compound_t *vc){

  unsigned int s;
  vrna_sc_t    *sc;

  if(vc){
    vrna_sc_remove(vc);

    switch(vc->type){
      case VRNA_VC_TYPE_SINGLE:     sc                    = (vrna_sc_t *)vrna_alloc(sizeof(vrna_sc_t));
                                    sc->energy_up         = NULL;
                                    sc->energy_bp         = NULL;
                                    sc->energy_stack      = NULL;
                                    sc->exp_energy_stack  = NULL;
                                    sc->exp_energy_up     = NULL;
                                    sc->exp_energy_bp     = NULL;
                                    sc->f                 = NULL;
                                    sc->exp_f             = NULL;
                                    sc->data              = NULL;
                                    sc->free_data         = NULL;

                                    vc->sc  = sc;
                                    break;

      case VRNA_VC_TYPE_ALIGNMENT:  vc->scs = (vrna_sc_t **)vrna_alloc(sizeof(vrna_sc_t*) * (vc->n_seq + 1));
                                    for(s = 0; s < vc->n_seq; s++){
                                      sc                    = (vrna_sc_t *)vrna_alloc(sizeof(vrna_sc_t));
                                      sc->energy_up         = NULL;
                                      sc->energy_bp         = NULL;
                                      sc->energy_stack      = NULL;
                                      sc->exp_energy_stack  = NULL;
                                      sc->exp_energy_up     = NULL;
                                      sc->exp_energy_bp     = NULL;
                                      sc->f                 = NULL;
                                      sc->exp_f             = NULL;
                                      sc->data              = NULL;
                                      sc->free_data         = NULL;

                                      vc->scs[s]  = sc;
                                    }
                                    break;
      default:                      /* do nothing */
                                    break;
    }
  }
}

PUBLIC void
vrna_sc_remove(vrna_fold_compound_t *vc){

  int s;

  if(vc){
    switch(vc->type){
      case  VRNA_VC_TYPE_SINGLE:    vrna_sc_free(vc->sc);
                                    vc->sc = NULL;
                                    break;
      case  VRNA_VC_TYPE_ALIGNMENT: if(vc->scs){
                                      for(s = 0; s < vc->n_seq; s++)
                                        vrna_sc_free(vc->scs[s]);
                                      free(vc->scs);
                                    }
                                    vc->scs = NULL;
                                    break;
      default:                      /* do nothing */
                                    break;
    }
  }
}

PUBLIC void
vrna_sc_free(vrna_sc_t *sc){

  int i;
  if(sc){
    if(sc->energy_up){
      for(i = 0; sc->energy_up[i]; free(sc->energy_up[i++]));
      free(sc->energy_up);
    }
    if(sc->exp_energy_up){
      for(i = 0; sc->exp_energy_up[i]; free(sc->exp_energy_up[i++]));
      free(sc->exp_energy_up);
    }
    
    free(sc->energy_bp);
    free(sc->exp_energy_bp);
    free(sc->energy_stack);
    free(sc->exp_energy_stack);

    if(sc->free_data)
      sc->free_data(sc->data);

    free(sc);
  }
}

PUBLIC void
vrna_sc_add_bp(vrna_fold_compound_t *vc,
                        const FLT_OR_DBL **constraints,
                        unsigned int options){
                        

  if(options & VRNA_CONSTRAINT_SOFT_MFE)
    sc_add_bp_mfe(vc, constraints, options);

  if(options & VRNA_CONSTRAINT_SOFT_PF)
    sc_add_bp_pf(vc, constraints, options);
}


PUBLIC  int
vrna_sc_add_SHAPE_zarringhalam( vrna_fold_compound_t *vc,
                                const double *reactivities,
                                double b,
                                double default_value,
                                const char *shape_conversion,
                                unsigned int options){

  int       i, j, n, ret;
  double        *pr;
  FLT_OR_DBL    *up, **bp;
  vrna_md_t *md;

  ret = 0; /* error */

  if(vc && reactivities && (vc->type == VRNA_VC_TYPE_SINGLE)){
    n   = vc->length;
    md  = (options & VRNA_CONSTRAINT_SOFT_PF) ? &(vc->exp_params->model_details) : &(vc->params->model_details);

    /* first we copy over the reactivities to convert them into probabilities later on */
    pr = (double *)vrna_alloc(sizeof(double) * (n + 1));
    for(i=0; i<=n; i++)
      pr[i] = reactivities[i];

    if(vrna_sc_SHAPE_to_pr(shape_conversion, pr, n, default_value)){

      /*  now, convert them into pseudo free energies for unpaired, and
          paired nucleotides
      */
      up = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
      bp = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 1));
      for(i = 1; i <= n; ++i){
        up[i] = b * fabs(pr[i] - 1);
        bp[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 1));
        for(j = i + md->min_loop_size + 1; j <= n; ++j)
          bp[i][j] = b * (pr[i] + pr[j]);
      }

      /* add the pseudo energies as soft constraints */
      vrna_sc_add_up(vc, (const FLT_OR_DBL *)up, options);
      vrna_sc_add_bp(vc, (const FLT_OR_DBL **)bp, options);

      /* clean up memory */
      for(i = 1; i <= n; ++i)
        free(bp[i]);
      free(bp);
      free(up);

      ret = 1; /* success */
    }

    free(pr);

  }

  return ret;
}


PUBLIC int
vrna_sc_add_SHAPE_deigan( vrna_fold_compound_t *vc,
                          const double *reactivities,
                          double m,
                          double b,
                          unsigned int options){

  int     i;
  FLT_OR_DBL  *values;

  if(vc && reactivities && (vc->type == VRNA_VC_TYPE_SINGLE)){

    values = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));

    /* first convert the values according to provided slope and intercept values */
    for (i = 1; i <= vc->length; ++i){
      values[i] = reactivities[i] < 0 ? 0. : (FLT_OR_DBL)(m * log(reactivities[i] + 1) + b);
    }

    if(options & VRNA_CONSTRAINT_SOFT_MFE)
      sc_add_stack_en_mfe(vc, (const FLT_OR_DBL *)values, options);

    if(options & VRNA_CONSTRAINT_SOFT_PF)
      sc_add_stack_en_pf(vc, (const FLT_OR_DBL *)values, options);

    free(values);
    return 1; /* success */
  } else {
    return 0; /* error */
  }
}

PUBLIC int
vrna_sc_add_SHAPE_deigan_ali( vrna_fold_compound_t *vc,
                              const char **shape_files,
                              const int *shape_file_association,
                              double m,
                              double b,
                              unsigned int options){

  float           reactivity, *reactivities, e1;
  char            *line, nucleotide, *sequence;
  int             s, i, p, r, position, *pseudo_energies, n_seq;
  unsigned short  **a2s;

  if(vc->type == VRNA_VC_TYPE_ALIGNMENT){
    n_seq = vc->n_seq;
    a2s   = vc->a2s;

    vrna_sc_init(vc);

    for(s = 0; shape_file_association[s] != -1; s++){
      int ss = shape_file_association[s]; /* actual sequence number in alignment */

      if(ss >= n_seq){
        vrna_message_warning("SHAPE file association exceeds sequence number in alignment");
        continue;
      }

      /* read the shape file */
      FILE *fp;
      if(!(fp = fopen(shape_files[s], "r"))){
        fprintf(stderr, "WARNING: SHAPE data file %d could not be opened. No shape data will be used.\n", s);
      } else {

        reactivities  = (float *)vrna_alloc(sizeof(float) * (vc->length + 1));
        sequence      = (char *)vrna_alloc(sizeof(char) * (vc->length + 1));

        /* initialize reactivities with missing data for entire alignment length */
        for(i = 1; i <= vc->length; i++)
          reactivities[i] = -1.;

        while((line=get_line(fp))){
          r = sscanf(line, "%d %c %f", &position, &nucleotide, &reactivity);
          if(r){
            if((position <= 0) || (position > vc->length))
              vrna_message_error("provided shape data outside of sequence scope");

            switch(r){
              case 1:   nucleotide = 'N';
                        /* fall through */
              case 2:   reactivity = -1.;
                        /* fall through */
              default:  sequence[position-1]    = nucleotide;
                        reactivities[position]  = reactivity;
                        break;
            }
          }
          free(line);
        }
        fclose(fp);

        sequence[vc->length] = '\0';

        /* double check information by comparing the sequence read from */
        char *tmp_seq = get_ungapped_sequence(vc->sequences[shape_file_association[s]]);
        if(strcmp(tmp_seq, sequence)){
          fprintf(stderr, "WARNING: Input sequence %d differs from sequence provided via SHAPE file!\n", shape_file_association[s]);
        }
        free(tmp_seq);

        /* convert reactivities to pseudo energies */
        for(i = 1; i <= vc->length; i++){
          if(reactivities[i] < 0)
            reactivities[i] = 0.;
          else
            reactivities[i] = m * log(reactivities[i] + 1.) + b; /* this should be a value in kcal/mol */
        }

        /*  begin actual storage of the pseudo energies */
        /*  beware of the fact that energy_stack will be accessed through a2s[s] array,
            hence pseudo_energy might be gap-free (default)
        */
        if(options & VRNA_CONSTRAINT_SOFT_MFE){
          int energy, cnt, gaps, is_gap;
          pseudo_energies = (int *)vrna_alloc(sizeof(int) * (vc->length + 1));
          for(gaps = cnt = 0, i = 1; i<=vc->length; i++){
            is_gap  = (vc->sequences[ss][i-1] == '-') ? 1 : 0;
            energy  = ((i - gaps > 0) && !(is_gap)) ? (int)(reactivities[i - gaps] * 100.) : 0;

            if(vc->params->model_details.oldAliEn){
              pseudo_energies[i] = energy;
              cnt++;
            } else if(!is_gap){ /* store gap-free */
              pseudo_energies[a2s[ss][i]] = energy;
              cnt++;
            }

            gaps += is_gap;
          }

          /* resize to actual number of entries */
          pseudo_energies = vrna_realloc(pseudo_energies, sizeof(int) * (cnt + 2));
          vc->scs[ss]->energy_stack = pseudo_energies;
        }

        if(options & VRNA_CONSTRAINT_SOFT_PF){
          FLT_OR_DBL *exp_pe = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));
          for(i=0;i<=vc->length;i++)
            exp_pe[i] = 1.;

          for(p = 0, i = 1; i<=vc->length; i++){
            e1 = (i - p > 0) ? reactivities[i - p] : 0.;
            if(vc->sequences[ss][i-1] == '-'){
              p++; e1 = 0.;
            }
            exp_pe[i] = (FLT_OR_DBL)exp(-(e1 * 1000.) / vc->exp_params->kT );
          }
          vc->scs[ss]->exp_energy_stack = exp_pe;
        }
        
        free(reactivities);
      }
    }

    return 1; /* success */
  } else {
    return 0; /* error */
  }
}

PUBLIC  int
vrna_sc_SHAPE_parse_method( const char *method_string,
                            char *method,
                            float *param_1,
                            float *param_2){

  const char *params = method_string + 1;

  *param_1 = 0;
  *param_2 = 0;

  if (!method_string || !method_string[0])
    return 0;

  *method = method_string[0];

  switch(method_string[0]){
    case 'Z':   *param_1 = 0.89;
                sc_parse_parameters(params, 'b', '\0', param_1, NULL);
                break;

    case 'D':   *param_1 = 1.8;
                *param_2 = -0.6;
                sc_parse_parameters(params, 'm', 'b', param_1, param_2);
                break;

    case 'W':   break;

    default:    *method = 0;
                return 0;
  }

  return 1;
}

PUBLIC void
vrna_sc_add_up(vrna_fold_compound_t *vc,
                        const FLT_OR_DBL *constraints,
                        unsigned int options){

  if(options & VRNA_CONSTRAINT_SOFT_MFE)
    sc_add_up_mfe(vc, constraints, options);

  if(options & VRNA_CONSTRAINT_SOFT_PF)
    sc_add_up_pf(vc, constraints, options);
}

PUBLIC void
vrna_sc_add_data( vrna_fold_compound_t *vc,
                  void *data,
                  vrna_callback_free_auxdata *free_data){

  if(vc){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->data        = data;
      vc->sc->free_data   = free_data;
    }
  }
}

PUBLIC void
vrna_sc_add_f(vrna_fold_compound_t *vc,
              vrna_callback_sc_energy *f){

  if(vc && f){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->f       = f;
    }
  }
}

PUBLIC void
vrna_sc_add_bt( vrna_fold_compound_t *vc,
                vrna_callback_sc_backtrack *f){

  if(vc && f){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->bt      = f;
    }
  }
}

PUBLIC void
vrna_sc_add_exp_f(vrna_fold_compound_t *vc,
                  vrna_callback_sc_exp_energy *exp_f){

  if(vc && exp_f){
    if(vc->type == VRNA_VC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->exp_f   = exp_f;
    }
  }
}

PRIVATE void
sc_parse_parameters( const char *string,
                        char c1,
                        char c2,
                        float *v1,
                        float *v2){

  char fmt[8];
  const char warning[] = "SHAPE method parameters not recognized! Using default parameters!";
  int r;

  assert(c1);
  assert(v1);

  if(!string || !(*string))
    return;

  if(c2 == 0 || v2 == NULL){
    sprintf(fmt, "%c%%f", c1);
    r = sscanf(string, fmt, v1);

    if(!r)
      vrna_message_warning(warning);

    return;
  }

  sprintf(fmt, "%c%%f%c%%f", c1, c2);
  r = sscanf(string, fmt, v1, v2);

  if(r!=2){
    sprintf(fmt, "%c%%f", c1);
    r = sscanf(string, fmt, v1);

    if(!r){
      sprintf(fmt, "%c%%f", c2);
      r = sscanf(string, fmt, v2);

      if(!r)
        vrna_message_warning(warning);
    }
  }
}

PRIVATE void
sc_add_bp_mfe(vrna_fold_compound_t *vc,
              const FLT_OR_DBL **constraints,
              unsigned int options){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;
  int           *idx;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_init(vc);

    sc              = vc->sc;
    sc->energy_bp = (int *)vrna_alloc(sizeof(int) * (((n + 1) * (n + 2)) / 2));

    idx = vc->jindx;
    for(i = 1; i < n; i++)
      for(j=i+1; j<=n; j++)
        sc->energy_bp[idx[j]+i] = (int)(constraints[i][j] * 100.);

  }
}

PRIVATE void
sc_add_bp_pf( vrna_fold_compound_t *vc,
              const FLT_OR_DBL **constraints,
              unsigned int options){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;
  int           *idx;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_init(vc);

    sc = vc->sc;

    vrna_exp_param_t  *exp_params = vc->exp_params;
    double            GT          = 0.;
    double            temperature = exp_params->temperature;
    double            kT          = exp_params->kT;
    double            TT          = (temperature+K0)/(Tmeasure);

    if(sc->exp_energy_bp)
      free(sc->exp_energy_bp);
    sc->exp_energy_bp     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (((n + 1) * (n + 2)) / 2));

    idx = vc->iindx;
    for(i = 1; i < n; i++)
      for(j=i+1; j<=n; j++){
        GT = constraints[i][j] * TT * 1000.;
        sc->exp_energy_bp[idx[i]-j] = (FLT_OR_DBL)exp( -GT / kT);
      }
  }
}

PRIVATE void
sc_add_up_mfe(vrna_fold_compound_t *vc,
              const FLT_OR_DBL *constraints,
              unsigned int options){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_init(vc);

    sc  = vc->sc;
    /*  allocate memory such that we can access the soft constraint
        energies of a subsequence of length j starting at position i
        via sc->energy_up[i][j]
    */
    if(sc->energy_up){
      for(i = 0; i <= n; i++)
        if(sc->energy_up[i])
          free(sc->energy_up[i]);
      free(sc->energy_up);
    }

    sc->energy_up = (int **)vrna_alloc(sizeof(int *) * (n + 2));
    for(i = 0; i <= n; i++)
      sc->energy_up[i] = (int *)vrna_alloc(sizeof(int) * (n - i + 2));

    sc->energy_up[n+1] = NULL;

    for(i = 1; i <= n; i++){
      for(j = 1; j <= (n - i + 1); j++){
        sc->energy_up[i][j] =   sc->energy_up[i][j-1]
                                  + (int)(constraints[i+j-1] * 100); /* convert to 10kal/mol */
      }
    }
  }
}

PRIVATE void
sc_add_up_pf( vrna_fold_compound_t *vc,
              const FLT_OR_DBL *constraints,
              unsigned int options){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_init(vc);

    sc  = vc->sc;

    vrna_exp_param_t   *exp_params = vc->exp_params;
    double             GT          = 0.;
    double             temperature = exp_params->temperature;
    double             kT          = exp_params->kT;
    double             TT          = (temperature+K0)/(Tmeasure);

    /* #################################### */
    /* # single nucleotide contributions  # */
    /* #################################### */

    /*  allocate memory such that we can access the soft constraint
        energies of a subsequence of length j starting at position i
        via sc->exp_energy_up[i][j]
    */
    if(sc->exp_energy_up){
      for(i = 0; i <= n; i++)
        if(sc->exp_energy_up[i])
          free(sc->exp_energy_up[i]);
      free(sc->exp_energy_up);
    }

    sc->exp_energy_up = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 2));
    for(i = 0; i <= n; i++){
      sc->exp_energy_up[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n - i + 2));
      for(j = 0; j < n - i + 2; j++)
        sc->exp_energy_up[i][j] = 1.;
    }

    sc->exp_energy_up[n+1] = NULL;

    for(i = 1; i <= n; i++){
      for(j = 1; j <= (n - i + 1); j++){
        GT  = (double)((int)(constraints[i+j-1] * 100)) * TT * 10.; /* convert to cal/mol */
        sc->exp_energy_up[i][j] =   sc->exp_energy_up[i][j-1]
                                      * (FLT_OR_DBL)exp( -GT / kT);
      }
    }
  }
}

PRIVATE void
sc_add_stack_en_mfe(vrna_fold_compound_t *vc,
                    const FLT_OR_DBL *constraints,
                    unsigned int options){
  int i;

  if(!vc->sc)
    vrna_sc_init(vc);

  if(!vc->sc->energy_stack)
    vc->sc->energy_stack = (int *)vrna_alloc(sizeof(int) * (vc->length + 1));

  for(i = 1; i <= vc->length; ++i)
    vc->sc->energy_stack[i] += (int)(constraints[i] * 100.);
}

PRIVATE void
sc_add_stack_en_pf( vrna_fold_compound_t *vc,
                    const FLT_OR_DBL *constraints,
                    unsigned int options){
  int i;

  if(!vc->sc)
    vrna_sc_init(vc);

  if(!vc->sc->exp_energy_stack){
    vc->sc->exp_energy_stack = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (vc->length + 1));
    for(i = 0; i <= vc->length; ++i)
      vc->sc->exp_energy_stack[i] = 1.;
  }

  for(i = 1; i <= vc->length; ++i)
    vc->sc->exp_energy_stack[i] *= (FLT_OR_DBL)exp(-(constraints[i] * 1000.)/ vc->exp_params->kT);
}

PRIVATE INLINE  void
adjust_ptypes(char *ptype,
              vrna_hc_t *hc,
              unsigned int length,
              unsigned int idx_type){

  unsigned  int i,j;
  int           *index;
  char          *matrix;

  matrix = hc->matrix;

  if(idx_type){
    index = vrna_idx_row_wise(length);
    for(i = 1; i < length; i++)
      for(j = i + 1; j <= length; j++)
        if(matrix[index[i] - j])
          if(!ptype[index[i] - j])
            ptype[index[i] - j] = 7; /* set to non-canonical pair */

  } else {
    index = vrna_idx_col_wise(length);
    for(i = 1; i < length; i++)
      for(j = i + 1; j <= length; j++)
        if(matrix[index[j] + i])
          if(!ptype[index[j] + i])
            ptype[index[j] + i] = 7; /* set to non-canonical pair */

  }
  free(index);
}
  
#ifdef  VRNA_BACKWARD_COMPAT

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

PUBLIC  void
print_tty_constraint_full(void){

  vrna_message_constraint_options_all();
}

PUBLIC  void
print_tty_constraint(unsigned int option){

  vrna_message_constraint_options(option);
}

PUBLIC void
constrain_ptypes( const char *constraint,
                  unsigned int length,
                  char *ptype,
                  int *BP,
                  int min_loop_size,
                  unsigned int idx_type){

  int n,i,j,k,l;
  int hx, *stack;
  char type;
  int *index;

  if(constraint == NULL) return;

  n = (int)strlen(constraint);

  stack = vrna_alloc(sizeof(int)*(n+1));

  if(!idx_type){ /* index allows access in energy matrices at pos (i,j) via index[j]+i */
    index = vrna_idx_col_wise(length);

    for(hx=0, j=1; j<=n; j++){
      switch(constraint[j-1]){
        case '|':   if(BP) BP[j] = -1;
                    break;
        case 'x':   /* can't pair */
                    for (l=1; l<j-min_loop_size; l++)
                      ptype[index[j]+l] = 0;
                    for (l=j+min_loop_size+1; l<=(int)length; l++)
                      ptype[index[l]+j] = 0;
                    break;
        case '(':   stack[hx++]=j;
                    /* fallthrough */
        case '<':   /* pairs upstream */
                    for (l=1; l<j-min_loop_size; l++)
                      ptype[index[j]+l] = 0;
                    break;
        case ')':   if (hx<=0) {
                      fprintf(stderr, "%s\n", constraint);
                      vrna_message_error("unbalanced brackets in constraint");
                    }
                    i = stack[--hx];
                    type = ptype[index[j]+i];
                    for (k=i+1; k<=(int)length; k++)
                      ptype[index[k]+i] = 0;
                    /* don't allow pairs i<k<j<l */
                    for (l=j; l<=(int)length; l++)
                      for (k=i+1; k<=j; k++)
                        ptype[index[l]+k] = 0;
                    /* don't allow pairs k<i<l<j */
                    for (l=i; l<=j; l++)
                      for (k=1; k<=i; k++)
                        ptype[index[l]+k] = 0;
                    for (k=1; k<j; k++)
                      ptype[index[j]+k] = 0;
                    ptype[index[j]+i] = (type==0) ? 7 : type;
                    /* fallthrough */
        case '>':   /* pairs downstream */
                    for (l=j+min_loop_size+1; l<=(int)length; l++)
                      ptype[index[l]+j] = 0;
                    break;
      }
    }
  }
  else{ /* index allows access in energy matrices at pos (i,j) via index[i]-j */
    index = vrna_idx_row_wise(length);

    for(hx=0, j=1; j<=n; j++) {
      switch (constraint[j-1]) {
        case 'x':   /* can't pair */
                    for (l=1; l<j-min_loop_size; l++)
                      ptype[index[l]-j] = 0;
                    for (l=j+min_loop_size+1; l<=(int)length; l++)
                      ptype[index[j]-l] = 0;
                    break;
        case '(':   stack[hx++]=j;
                    /* fallthrough */
        case '<':   /* pairs upstream */
                    for (l=1; l<j-min_loop_size; l++)
                      ptype[index[l]-j] = 0;
                    break;
        case ')':   if (hx<=0) {
                      fprintf(stderr, "%s\n", constraint);
                      vrna_message_error("unbalanced brackets in constraints");
                    }
                    i = stack[--hx];
                    type = ptype[index[i]-j];
                    /* don't allow pairs i<k<j<l */
                    for (k=i; k<=j; k++)
                      for (l=j; l<=(int)length; l++)
                        ptype[index[k]-l] = 0;
                    /* don't allow pairs k<i<l<j */
                    for (k=1; k<=i; k++)
                      for (l=i; l<=j; l++)
                        ptype[index[k]-l] = 0;
                    ptype[index[i]-j] = (type==0) ? 7 : type;
                    /* fallthrough */
        case '>':   /* pairs downstream */
                    for (l=j+min_loop_size+1; l<=(int)length; l++)
                      ptype[index[j]-l] = 0;
                    break;
      }
    }
  }
  if (hx!=0) {
    fprintf(stderr, "%s\n", constraint);
    vrna_message_error("unbalanced brackets in constraint string");
  }
  free(index);
  free(stack);
}

#endif
