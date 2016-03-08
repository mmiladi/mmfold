/**********************************************/
/* BEGIN interface for Suboptimal Structure   */
/* predictionn                                */
/**********************************************/

// from subopt.h


extern  int subopt_sorted;                       /* sort output by energy */

typedef struct {
  float energy;
  char  *structure;
} SOLUTION;

%extend SOLUTION {
        SOLUTION *get(int i) {
//           static int size=-1;
//           if (size<0) {
//             SOLUTION *s;
//             for (s=self; s->structure; s++);
//             size= (int) (s-self);
//           }
//           if (i>=size) {
//             warn("value out of range");
//             return NULL;
//           }
           return self+i;
        }

        int size() {
           SOLUTION *s;
           for (s=self; s->structure; s++);
           return (int)(s-self);
        }

        ~SOLUTION() {
           SOLUTION *s;
           for (s=self; s->structure; s++) free(s->structure);
           free(self);
        }
}

%newobject subopt;
extern  SOLUTION *subopt (char *seq, char *constraint, int delta, FILE *fp=NULL);

%ignore subopt_par;
%ignore subopt_circ;
%ignore zukersubopt;
%ignore zukersubopt_par;

%include  <ViennaRNA/subopt.h>
