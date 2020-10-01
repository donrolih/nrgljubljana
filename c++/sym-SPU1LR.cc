class SymmetrySPU1LR : public SymFieldLR {
  private:
  outfield Sz2, Sz;

  public:
  SymmetrySPU1LR() : SymFieldLR() { all_syms["SPU1LR"] = this; }

  void init() override {
    Sz2.set("<Sz^2>", 1);
    Sz.set("<Sz>", 2);
    InvarStructure InvStruc[] = {
       {"SSZ", additive},    // spin projection
       {"P", multiplicative} // parity
    };
    initInvar(InvStruc, ARRAYLENGTH(InvStruc));
    InvarSinglet = Invar(0, 1);
  }

  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) override {
    // The spin projection of the operator is defined by the difference
    // in Sz of both the invariant subspaces.
    SZspin ssz1  = I1.get("SSZ");
    SZspin sszp  = Ip.get("SSZ");
    SZspin sszop = ssz1 - sszp;
    return sszop == SPIN;
  }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) override {
    return u1_equality(I1.get("SSZ"), I2.get("SSZ"), I3.get("SSZ")) && z2_equality(I1.get("P"), I2.get("P"), I3.get("P"));
  }

  void load() override {
    switch (channels) {
      case 1:
#include "spu1lr/spu1lr-1ch-In2.dat"
#include "spu1lr/spu1lr-1ch-QN.dat"
        break;

      case 2:
#include "spu1lr/spu1lr-2ch-In2.dat"
#include "spu1lr/spu1lr-2ch-QN.dat"
        break;

      default: my_assert_not_reached();
    }
  }

  void calculate_TD(const DiagInfo &diag, double factor) override {
    bucket trSZ, trSZ2; // Tr[S_z], Tr[S_z^2]

    for (const auto &[I, eig]: diag) {
      const SZspin ssz  = I.get("SSZ");
      const double sumZ = calculate_Z(I, eig, factor);

      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * sqr(SZ(ssz));
    }

    Sz2 = trSZ2 / stats.Z;
    Sz  = trSZ / stats.Z;
  }

  DECL;
  HAS_DOUBLET;
  HAS_TRIPLET;
  HAS_GLOBAL;
};

Symmetry *SymSPU1LR = new SymmetrySPU1LR;

#undef ISOSPINX
#define ISOSPINX(i, j, ch, factor) diag_offdiag_function(i, j, ch, t_matel(factor) * 2.0 * delta(step.N() + 1, ch), h, qq)

#undef ANOMALOUS
#define ANOMALOUS(i, j, ch, factor) offdiag_function(i, j, ch, 0, t_matel(factor) * kappa(step.N(), ch), h, qq, In, opch)

#undef OFFDIAG
#define OFFDIAG(i, j, ch, factor0) offdiag_function(i, j, ch, 0, t_matel(factor0) * xi(step.N(), ch), h, qq, In, opch)

#undef DIAG
#define DIAG(i, ch, number) diag_function(i, ch, number, zeta(step.N() + 1, ch), h, qq)

void SymmetrySPU1LR::makematrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, const Opch &opch) {
  switch (channels) {
    case 1:
#include "spu1lr/spu1lr-1ch-offdiag.dat"
#include "spu1lr/spu1lr-1ch-anomalous.dat"
#include "spu1lr/spu1lr-1ch-diag.dat"
#include "spu1lr/spu1lr-1ch-isospinx.dat"
      break;

    case 2:
#include "spu1lr/spu1lr-2ch-diag.dat"
#include "spu1lr/spu1lr-2ch-offdiag.dat"
#include "spu1lr/spu1lr-2ch-anomalous.dat"
#include "spu1lr/spu1lr-2ch-isospinx.dat"
      break;

    default: my_assert_not_reached();
  }
}

#include "nrg-recalc-SPU1LR.cc"
