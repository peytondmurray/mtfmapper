#ifndef LOESS_PARMS_H
#define LOESS_PARMS_H

class Loess_parms {
  public:
      Loess_parms(double alpha=0.5) : alpha(alpha) {}
      
      void set_alpha(double a) {
          alpha = a;
      }
        
      double get_alpha(void) const {
          return alpha;
      }
      
      static Loess_parms& get_instance(void) {
          static Loess_parms singleton;
          return singleton;
      }
      
  private:
      double alpha = 0.5;
};


#endif
