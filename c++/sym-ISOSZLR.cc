template<typename SC>
class SymmetryISOSZLR_tmpl : public SymFieldLR_tmpl<SC> {
 private:
   outfield Sz2, Sz, Q2;
   using Symmetry_tmpl<SC>::P;
   using Symmetry_tmpl<SC>::In;
   using Symmetry_tmpl<SC>::QN;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryISOSZLR_tmpl(const Params &P, Allfields &allfields) : SymFieldLR_tmpl<SC>(P),
     Sz2(P, allfields, "<Sz^2>", 1), Sz(P, allfields, "<Sz>", 2), Q2(P, allfields, "<Q^2>", 3) {
       initInvar({
         {"II", additive},     // isospin
         {"SSZ", additive},    // spin projection
         {"P", multiplicative} // parity
       });
       this->InvarSinglet = Invar(1, 0, 1);
     }

  // Multiplicity of the I=(II,SSZ) subspace = (2I+1) = II.
  size_t mult(const Invar &I) const override { return I.get("II"); }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    SZspin ssz1  = I1.get("SSZ");
    SZspin sszp  = Ip.get("SSZ");
    SZspin sszop = ssz1 - sszp;
    return sszop == SPIN;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
    return u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")) && su2_triangle_inequality(I1.get("II"), I2.get("II"), I3.get("II"))
       && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    my_assert(P.channels == 2);
#include "isoszlr/isoszlr-2ch-In2.dat"
#include "isoszlr/isoszlr-2ch-QN.dat"
  }

  double specdens_factor(const Invar &Ip, const Invar &I1) const override {
    check_abs_diff(Ip, I1, "SSZ", 1);
    const Ispin iip = Ip.get("II");
    const Ispin ii1 = I1.get("II");
    const double isofactor = (ii1 == iip + 1 ? ISO(iip) + 1.0 : ISO(iip));
    return isofactor;
  }

  void calculate_TD(const Step &step, const DiagInfo_tmpl<SC> &diag, const Stats_tmpl<SC> &stats, const double factor) override {
    bucket trSZ, trSZ2, trIZ2; // Tr[S_z], Tr[S_z^2], Tr[I_z^2]
    for (const auto &[I, eig]: diag) {
      const Ispin ii    = I.get("II");
      const SZspin ssz  = I.get("SSZ");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * pow(SZ(ssz),2);        // isospin multiplicity contained in sumZ
      trIZ2 += sumZ * (ii * ii - 1) / 12.; // spin multiplicity contained in sumZ
    }
    Sz  = trSZ / stats.Z;
    Sz2 = trSZ2 / stats.Z;
    Q2  = (4 * trIZ2) / stats.Z;
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
};

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(step, i, j, ch, 0, t_matel(factor0) * coef.xi(step.N(), ch), h, qq, In, opch)

template<typename SC>
void SymmetryISOSZLR_tmpl<SC>::make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch_tmpl<SC> &opch, const Coef_tmpl<SC> &coef) {
  Ispin ii = I.get("II");
#include "isoszlr/isoszlr-2ch-offdiag.dat"
}

#include "nrg-recalc-ISOSZLR.cc"
