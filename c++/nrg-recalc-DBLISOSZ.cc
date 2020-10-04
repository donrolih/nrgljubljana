// *** WARNING!!! Modify nrg-recalc-DBLISOSZ.cc.m4, not nrg-recalc-DBLISOSZ.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Mar 2010
// This file pertains to (I1,I2,Sz) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers







  





namespace DBLISOSZ {
#include "dblisosz/dblisosz-2ch-def.dat"
}

// Recalculate matrix elements of a doublet tenzor operator
MatrixElements SymmetryDBLISOSZ::recalc_doublet(const DiagInfo &diag, const QSrmax &qsrmax, const MatrixElements &cold) {
  MatrixElements cnew;
  for(const auto &[I1, eig]: diag) {
    Ispin ii11 = I1.get("II1");
    Ispin ii21 = I1.get("II2");
    Sspin ssz1 = I1.get("SSZ");
    Invar Ip;

    Ip = Invar(ii11 - 1, ii21, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletm0m.dat" << ", len=" << DBLISOSZ::LENGTH_D1_2CH << ", Iop=" << Invar(2, 1, +1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "dblisosz/dblisosz-2ch-doubletm0m.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_D1_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, DBLISOSZ::LENGTH_D1_2CH, Invar(2, 1, +1));
  }
};

    Ip = Invar(ii11 - 1, ii21, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletm0p.dat" << ", len=" << DBLISOSZ::LENGTH_D1_2CH << ", Iop=" << Invar(2, 1, -1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "dblisosz/dblisosz-2ch-doubletm0p.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_D1_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, DBLISOSZ::LENGTH_D1_2CH, Invar(2, 1, -1));
  }
};

    Ip = Invar(ii11 + 1, ii21, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletp0m.dat" << ", len=" << DBLISOSZ::LENGTH_D1_2CH << ", Iop=" << Invar(2, 1, +1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "dblisosz/dblisosz-2ch-doubletp0m.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_D1_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, DBLISOSZ::LENGTH_D1_2CH, Invar(2, 1, +1));
  }
};

    Ip = Invar(ii11 + 1, ii21, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doubletp0p.dat" << ", len=" << DBLISOSZ::LENGTH_D1_2CH << ", Iop=" << Invar(2, 1, -1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "dblisosz/dblisosz-2ch-doubletp0p.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_D1_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, DBLISOSZ::LENGTH_D1_2CH, Invar(2, 1, -1));
  }
};

    Ip = Invar(ii11, ii21 - 1, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0mm.dat" << ", len=" << DBLISOSZ::LENGTH_D2_2CH << ", Iop=" << Invar(1, 2, +1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "dblisosz/dblisosz-2ch-doublet0mm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_D2_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, DBLISOSZ::LENGTH_D2_2CH, Invar(1, 2, +1));
  }
};

    Ip = Invar(ii11, ii21 - 1, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0mp.dat" << ", len=" << DBLISOSZ::LENGTH_D2_2CH << ", Iop=" << Invar(1, 2, -1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "dblisosz/dblisosz-2ch-doublet0mp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_D2_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, DBLISOSZ::LENGTH_D2_2CH, Invar(1, 2, -1));
  }
};

    Ip = Invar(ii11, ii21 + 1, ssz1 - 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0pm.dat" << ", len=" << DBLISOSZ::LENGTH_D2_2CH << ", Iop=" << Invar(1, 2, +1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "dblisosz/dblisosz-2ch-doublet0pm.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_D2_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, DBLISOSZ::LENGTH_D2_2CH, Invar(1, 2, +1));
  }
};

    Ip = Invar(ii11, ii21 + 1, ssz1 + 1);
    {
  nrglog('f', "RECALC(fn=" << "dblisosz/dblisosz-2ch-doublet0pp.dat" << ", len=" << DBLISOSZ::LENGTH_D2_2CH << ", Iop=" << Invar(1, 2, -1) << ")");
  if (diag.count(Ip)) {
    struct Recalc recalc_table[] = {
#include "dblisosz/dblisosz-2ch-doublet0pp.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_D2_2CH);
    recalc_general(diag, qsrmax, cold, cnew, I1, Ip, recalc_table, DBLISOSZ::LENGTH_D2_2CH, Invar(1, 2, -1));
  }
};
  }
  return cnew;
}

// Driver routine for recalc_f()
Opch SymmetryDBLISOSZ::recalc_irreduc(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, const Params &P) {
  Opch opch = newopch(P);
  for(const auto &[Ip, eig]: diag) {
    Invar I1;

    Ispin ii1p = Ip.get("II1");
    Ispin ii2p = Ip.get("II2");
    Sspin sszp = Ip.get("SSZ");

    // NN is index n of f_n, the last site in the chain prior to adding
    // the new site (f_{n+1}).
    int NN = step.getnn();

    // RECALC_F_TAB_... (filename, channel_number, array_length)

    I1 = Invar(ii1p + 1, ii2p, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isoup-a.dat" << ", ch=" << 0 << ", len=" << DBLISOSZ::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblisosz/dblisosz-2ch-type1-isoup-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, DBLISOSZ::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p + 1, ii2p, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isoup-a.dat" << ", ch=" << 0 << ", len=" << DBLISOSZ::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblisosz/dblisosz-2ch-type2-isoup-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, DBLISOSZ::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p, ii2p + 1, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isoup-b.dat" << ", ch=" << 1 << ", len=" << DBLISOSZ::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblisosz/dblisosz-2ch-type1-isoup-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[1][0], Ip, I1, recalc_table, DBLISOSZ::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p, ii2p + 1, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isoup-b.dat" << ", ch=" << 1 << ", len=" << DBLISOSZ::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblisosz/dblisosz-2ch-type2-isoup-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[1][0], Ip, I1, recalc_table, DBLISOSZ::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p - 1, ii2p, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isodown-a.dat" << ", ch=" << 0 << ", len=" << DBLISOSZ::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblisosz/dblisosz-2ch-type1-isodown-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, DBLISOSZ::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p - 1, ii2p, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isodown-a.dat" << ", ch=" << 0 << ", len=" << DBLISOSZ::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblisosz/dblisosz-2ch-type2-isodown-a.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[0][0], Ip, I1, recalc_table, DBLISOSZ::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p, ii2p - 1, sszp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type1-isodown-b.dat" << ", ch=" << 1 << ", len=" << DBLISOSZ::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblisosz/dblisosz-2ch-type1-isodown-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[1][0], Ip, I1, recalc_table, DBLISOSZ::LENGTH_I_2CH);
  }
};

    I1 = Invar(ii1p, ii2p - 1, sszp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "dblisosz/dblisosz-2ch-type2-isodown-b.dat" << ", ch=" << 1 << ", len=" << DBLISOSZ::LENGTH_I_2CH << ")");
  if (diag.count(I1)) {
    struct Recalc_f recalc_table[] = {
#include "dblisosz/dblisosz-2ch-type2-isodown-b.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == DBLISOSZ::LENGTH_I_2CH);
    recalc_f(diag, qsrmax, opch[1][0], Ip, I1, recalc_table, DBLISOSZ::LENGTH_I_2CH);
  }
};
  }
  return opch;
}

#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, qsrmax, I1, cn, i1, ip, value)

void SymmetryDBLISOSZ::recalc_global(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, string name, MatrixElements &cnew) {
  if (name == "SZtot") {
   for(const auto &[I1, eig]: diag) {
      const Twoinvar II{I1, I1};
      Matrix &cn        = cnew[II];
      switch (channels) {
        case 2:
#include "dblisosz/dblisosz-2ch-spinz.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}
