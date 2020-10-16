// symmetry.cc - Classes representing various symmetry types
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _symmetry_cc_
#define _symmetry_cc_

// Check if the triangle inequality is satisfied (i.e. if
// Clebsch-Gordan coefficient can be different from zero). This is
// important, for example, for triplet operators, which are zero when
// evaluated between two singlet states.
// Arguments ss1, ss2, ss3 are spin multiplicities.
// Returns true if the inequality is satisfied, false otherwise.
bool su2_triangle_inequality(int ss1, int ss2, int ss3) {
  return (abs(ss1 - ss2) <= ss3 - 1) && (abs(ss2 - ss3) <= ss1 - 1) && (abs(ss3 - ss1) <= ss2 - 1);
}

// The equality for U(1) symmetry.
bool u1_equality(int q1, int q2, int q3) { return q1 == q2 + q3; }

bool z2_equality(int p1, int p2, int p3) { return p1 == p2 * p3; }

// C_3 quantum number: Equality modulo 3
bool c3_equality(int p1, int p2, int p3) { return p1 == (p2 + p3) % 3; }

void opch1clear(Opch &opch, int i, const Params &P)
{
  opch[i].resize(P.perchannel);
  for (size_t j = 0; j < P.perchannel; j++) 
    opch[i][j].clear(); // set all ublas matrix elements to zero
}

Opch newopch(const Params &P)
{
  Opch opch(P.channels);
  for (size_t i = 0; i < P.channels; i++)
    opch1clear(opch, i, P);
  return opch;
}

struct Recalc_f {
    size_t i1; // subspace indexes
    size_t ip;
    t_factor factor;
};

// Structure which holds subspace information and factor for each of nonzero irreducible matrix elements. cf.
// Hofstetter PhD p. 120. <Q+1 S+-1/2 .. i1 ||f^\dag|| Q S .. ip>_N = factor < IN1 .. ||f^\dag|| INp ..>_{N_1}
struct Recalc {
    size_t i1{}; // combination of states
    size_t ip{};
    Invar IN1; // subspace in N-1 stage
    Invar INp;
    t_factor factor{}; // additional multiplicative factor
};

class Symmetry {
 protected:
   const Params &P;
   Allfields &allfields;

 public:
   // In and QN contain information about how the invariant subspaces at consecutive iteration steps are combined. In
   // is the array of quantum number DIFFERENCES used in the construction of the basis. QN is the array of the
   // conserved quantum numbers corresponding to the states being added. For example, in case of SU(2) symmetry, In
   // will include S_z, while QN will include S. In other words, QN involves those quantum numbers that we need to
   // retain in the calculation, while In involves those quantum numbers that "drop out" of the problem due to the
   // symmetry.
   std::vector<Invar> In, QN;

   Symmetry(const Params &P_, Allfields &allfields_) : P(P_), allfields(allfields_), In(P.combs+1), QN(P.combs+1) {}

   InvarVec input_subspaces() const { return In; }

   Invar ancestor(const Invar &I, const size_t i) const {
     const auto input = input_subspaces();
     Invar anc = I;
     anc.combine(input[i]);
     return anc; // I.combine(input[i]) == input[i].combine(I)
   }

   InvarVec ancestors(const Invar &I) const {
     auto input = input_subspaces();
     for (size_t i = 1; i <= P.combs; i++)
       input[i].combine(I); // In is the list of differences wrt I
     return input;
   }

   InvarVec dmnrg_subspaces(const Invar &I) const {
     auto input = input_subspaces();
     for (size_t i = 1; i <= P.combs; i++) {
       input[i].inverse();
       input[i].combine(I);
     }
     return input;
   }

   Invar InvarSinglet; // QNs for singlet operator
   Invar Invar_f;      // QNs for f operator

   size_t get_combs() const { return P.combs; }
   auto combs1() const { return boost::irange(size_t(1), size_t(P.combs+1)); } // note: 1-based!

   // For some symmetry types with two-channels we distinguish between
   // even and odd parity with respect to the channel-interchange
   // operation.
   virtual bool islr() { return false; }

   // Ditto for 3 channels: C_3 symmetry.
   virtual bool isc3() { return false; }

   // For some symmetry types, we may distinguish between spin-up
   // and spin-down quantities (in particular spin-up and spin-down
   // spectral functions).
   virtual bool isfield() { return false; }

   // Multiplicity of the states in the invariant subspace
   virtual size_t mult(const Invar &) const {
     return 1;
   };

   auto multfnc() const {
     return [this](const Invar &I) { return this->mult(I); };
   }

   double calculate_Z(const Invar &I, const Eigen &eig, double rescale_factor) const {
     return mult(I) * ranges::accumulate(eig.value_zero, 0.0, [rf=rescale_factor](auto sum, const auto &x) { return sum+exp(-rf*x); });
   }

   // Does the combination of invariant subspaces I1 and I2 contribute
   // to the spectral function corresponding to spin SPIN?
   virtual bool check_SPIN(const Invar &I1, const Invar &I2, const int &SPIN) const {
     my_assert(SPIN == 0);
     return true;
   }

   // Is the triangle inequality satisfied for I1, I2, and I3 (i.e. if
   // Clebsch-Gordan coefficient can be different from zero).  This is
   // important, for example, for triplet operators which are zero
   // when evaluated between two singlet states.
   virtual bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const { return true; }

   // Setup the combinations of quantum numbers that are used in the
   // construction of the Hamiltonian matrix.
   virtual void load() = 0;

   // Is an invariant subspace with given quantum numbers allowed?
   virtual bool Invar_allowed(const Invar &I) const { return true; }

   bool offdiag_contributes(const size_t i, const size_t j, const Rmaxvals &qq) const;
   void offdiag_function_impl(const Step &step, const size_t i, const size_t j, const size_t ch, const size_t fnr, const t_matel factor,
                              Matrix &h, const Rmaxvals &qq, const InvarVec &In, const Opch &opch) const;
   void diag_function_impl(const Step &step, const size_t i, const size_t ch, const double number, const t_coef sc_zeta, 
                           Matrix &h, const Rmaxvals &qq, const double f) const;
   void diag_function(const Step &step, const size_t i, const size_t ch, const double number, const t_coef sc_zeta, 
                      Matrix &h, const Rmaxvals &qq) const;
   void diag_function_half(const Step &step, const size_t i, const size_t ch, const double number, const t_matel sc_zeta,
                           Matrix &h, const Rmaxvals &qq) const;
   void diag_offdiag_function(const Step &step, const size_t i, const size_t j, const size_t chin, const t_matel factor,
                              Matrix &h, const Rmaxvals &qq) const;

   virtual void make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) = 0;

   // Called from recalc_dynamicsusceptibility().  This is the factor due
   // to the spin degeneracy when calculating the trace of Sz.Sz.
   virtual double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const { return 1.0; }

   // Called from recalc_dynamic_orb_susceptibility().  This is the factor due
   // to the orbital moment degeneracy when calculating the trace of Tz.Tz.
   virtual double dynamic_orb_susceptibility_factor(const Invar &Ip, const Invar &I1) const { return 1.0; }

   // Called from calc_specdens().
   // See spectral_density_clebschgordan.nb and DMNRG_clebschgordan.nb.
   virtual double specdens_factor(const Invar &Ip, const Invar &I1) const { return 1.0; }
   virtual double specdensquad_factor(const Invar &Ip, const Invar &I1) const { return 1.0; }

   virtual void calculate_TD(const Step &step, const DiagInfo &diag, const Stats &stats, double factor) = 0;

   virtual Opch recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) { my_assert_not_reached(); }
   virtual OpchChannel recalc_irreduc_substeps(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P, int M) { my_assert_not_reached(); }
   virtual MatrixElements recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) { my_assert_not_reached(); }
   virtual MatrixElements recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) { my_assert_not_reached(); }
   virtual MatrixElements recalc_orb_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) { my_assert_not_reached(); }
   virtual MatrixElements  recalc_quadruplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) { my_assert_not_reached(); }
   virtual void recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) { my_assert_not_reached(); }

   // Recalculates irreducible matrix elements of a singlet operator, as well as odd-parity spin-singlet operator (for
   //  parity -1). Generic implementation, valid for all symmetry types.
   MatrixElements recalc_singlet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &nold, int parity) {
     MatrixElements nnew;
     Recalc recalc_table[get_combs()];
     my_assert(islr() ? parity == 1 || parity == -1 : parity == 1);
     for (const auto &I : diag.subspaces()) {
       const Invar I1 = I;
       const Invar Ip = parity == -1 ? I.InvertParity() : I;
       for (size_t i = 1; i <= get_combs(); i++) {
         const auto anc = ancestor(I, i);
         recalc_table[i - 1] = {i, i, anc, parity == -1 ? anc.InvertParity() : anc, 1.0};
       }
       const auto Iop = parity == -1 ? InvarSinglet.InvertParity() : InvarSinglet;
       nnew[Twoinvar(I1,Ip)] = recalc_general(diag, qsrmax, nold, I1, Ip, recalc_table, get_combs(), Iop);
     }
     return nnew;
   }

   virtual void show_coefficients(const Step &step, const Coef &coef) {
     cout << setprecision(std::numeric_limits<double>::max_digits10);
     if (!P.substeps) {
       for (size_t i = 0; i < P.coefchannels; i++) {
         auto N = step.N();
         cout << "[" << i + 1 << "]"
           << " xi(" << N << ")=" << coef.xi(N, i) << " xi_scaled(" << N << ")=" << coef.xi(N, i)/step.scale()
             << " zeta(" << N+1 << ")=" << coef.zeta(N+1, i) << endl;
       }
     } else {
       const auto [N, M] = step.NM();
       for (auto i = 0; i < P.coeffactor; i++) {
         auto index = M + P.channels * i;
         cout << "[" << index << "]"
           << " xi(" << N << ")=" << coef.xi(N, index) << " zeta(" << N+1 << ")=" << coef.zeta(N+1, index) << endl;
       }
     }
   }

   virtual bool recalc_f_coupled(const Invar &I1, const Invar &I2, const Invar &If) { return true; } // used in recalc_f()

   Matrix recalc_f(const DiagInfo &diag, const QSrmax &qsrmax, const Invar &I1,
                   const Invar &Ip, const struct Recalc_f table[], const size_t jmax);

   Matrix recalc_general(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold,
                         const Invar &I1, const Invar &Ip, const struct Recalc table[], const size_t jmax, const Invar &Iop) const;

   void recalc1_global(const DiagInfo &diag, const QSrmax &qsrmax, const Invar &I,
                       Matrix &m, const size_t i1, const size_t ip, const t_factor value) const;

   auto CorrelatorFactorFnc() const   { return [this](const Invar &Ip, const Invar &I1) { return this->mult(I1); }; }
   auto SpecdensFactorFnc() const     { return [this](const Invar &Ip, const Invar &I1) { return this->specdens_factor(Ip, I1); }; }
   auto SpecdensquadFactorFnc() const { return [this](const Invar &Ip, const Invar &I1) { return this->specdensquad_factor(Ip, I1); }; }
   auto SpinSuscFactorFnc() const     { return [this](const Invar &Ip, const Invar &I1) { return this->dynamicsusceptibility_factor(Ip, I1); }; }
   auto OrbSuscFactorFnc() const      { return [this](const Invar &Ip, const Invar &I1) { return this->dynamic_orb_susceptibility_factor(Ip, I1); }; }
   auto TrivialCheckSpinFnc() const   { return [this](const Invar &Ip, const Invar &I1, int SPIN) { return true; }; }
   auto SpecdensCheckSpinFnc() const  { return [this](const Invar &I1, const Invar &Ip, int SPIN) { return this->check_SPIN(I1, Ip, SPIN); }; }
};

// Add DECL declaration in each symmetry class
#define DECL                                                                                                                                           \
  void make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch, const Coef &coef) override;  \
  Opch recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) override

// Optional declaration
#define HAS_SUBSTEPS OpchChannel recalc_irreduc_substeps(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P, int M) override
#define HAS_DOUBLET MatrixElements recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) override
#define HAS_TRIPLET MatrixElements recalc_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) override
#define HAS_ORB_TRIPLET MatrixElements recalc_orb_triplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) override
#define HAS_QUADRUPLET MatrixElements recalc_quadruplet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) override
#define HAS_GLOBAL void recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) override

class SymField : public Symmetry {
 public:
   template<typename ... Args> explicit SymField(Args&& ... args) : Symmetry(std::forward<Args>(args)...) {}
   bool isfield() override { return true; }
};

class SymLR : public Symmetry {
 public:
   template<typename ... Args> explicit SymLR(Args&& ... args) : Symmetry(std::forward<Args>(args)...) {}
   bool islr() override { return true; }
};

class SymC3 : public Symmetry {
 public:
   template<typename ... Args> explicit SymC3(Args&& ... args) : Symmetry(std::forward<Args>(args)...) {}
   bool isc3() override { return true; }
};

class SymFieldLR : public Symmetry {
  public:
   template<typename ... Args> explicit SymFieldLR(Args&& ... args) : Symmetry(std::forward<Args>(args)...) {}
   bool isfield() override { return true; }
   bool islr() override { return true; }
};

// Helper functions
void check_abs_diff(const Invar &Ip, const Invar &I1, const string &what, int diff) {
  const int a = Ip.get(what);
  const int b = I1.get(what);
  my_assert(abs(b - a) == diff);
}

void check_diff(const Invar &Ip, const Invar &I1, const string &what, int diff) {
  const int a = Ip.get(what);
  const int b = I1.get(what);
  my_assert(b - a == diff);
}

#endif // _symmetry_cc_
