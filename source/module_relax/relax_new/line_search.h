#ifndef LINE_SEARCH_H
#define LINE_SEARCH_H

class Line_Search
{
  public:
    Line_Search() {};
    ~Line_Search() {};

    bool line_search(const bool restart,     // whether to restart line search (when cg direction has changed)
                     const double x,         // position for function at x
                     const double y,         // value of function at x
                     const double f,         // gradient of function at x
                     double& xnew,           // postion where function has to be evaluated
                     const double conv_thr); // predicted change of function value of function at xnew

    bool first_order(const double x, const double y, const double f, double& xnew);

    bool third_order(const double x, const double y, const double f, double& xnew, const double conv_thr);

    void init_brent(const double x, const double y, const double f);

    void update_brent(const double x, const double y, const double f);

    bool brent(const double x, const double y, const double f, double& xnew, const double conv_thr);

  private:
    int ls_step = 0;      // step of line search
    bool bracked = false; // whether the minima is bracked by [a,c]
    /// these variables used to keep record of some points
    double xa = 0.0;
    double xb = 0.0;
    double xc = 0.0;
    double fa = 0.0;
    double fb = 0.0;
    double fc = 0.0;
    double ya = 0.0;
    double yb = 0.0;
    /// ------------------------------------
    double fstart = 0.0; // keep record of initial gradient for brent method
    const double e8 = 1.0e-8;
};

#endif