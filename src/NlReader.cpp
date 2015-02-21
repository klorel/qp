/* $Id: OSnl2OS.cpp 4210 2011-06-28 09:44:54Z stefan $ */
/** @file OSnl2OS.cpp
 * 
 * @author  Horand Gassmann, Jun Ma, Kipp Martin
 *
 * \remarks
 * Copyright (C) 2012, Horand Gassmann, Jun Ma, Kipp Martin,
 * Dalhousie University, and the University of Chicago.
 * All Rights Reserved.
 * This software is licensed under the Eclipse Public License. 
 * Please see the accompanying LICENSE file in root directory for terms.
 * 
 */

/**  
 *
 * <p><code>OSnl2OS</code> is used to convert information in AMPL nl format
 * to OS objects (OSInstance and OSOption/OSoL string) </p>
 *
 */

#include "NlReader.h"
#include "Problem.h"
#include <iostream>
std::map<size_t, std::string> NlReader::OppName =
    std::map<size_t, std::string>();
#include "nlp.h"
#include "getstub.h"
#include "r_opn.hd" /* for N_OPS */
#include "opcode.hd"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_CMATH_H
#  include <cmath.h>
# endif
#endif

#include <sstream>

#define R_OPS  ((ASL_fg*)asl)->I.r_ops_
#define OBJ_DE ((ASL_fg*)asl)->I.obj_de_
#define VAR_E  ((ASL_fg*)asl)->I.var_e_
#define CON_DE ((ASL_fg*)asl)->I.con_de_

#include <asl.h>
//#include "sufinfo.h"

using std::cerr;
using std::cout;
using std::endl;

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif

#define UN_IMPLEMENTTED_OPP(__OPP__) \
    std::cerr << "NOT IMPLEMENTED "<<#__OPP__ << __FILE__ << ",l:" << __LINE__ << std::endl
//#define AMPLDEBUG
void NlReader::free() {
  if (asl != 0) {
    std::free(A_vals);
    if (X0)
      delete[] X0;
    X0 = 0;
    if (havex0)
      delete[] havex0;
    havex0 = 0;
    if (pi0)
      delete[] pi0;
    pi0 = 0;
    if (havepi0)
      delete[] havepi0;
    havepi0 = 0;
    ASL_free(&asl);
  }
}
void NlReader::getCW() {
  if (!_isCW) {
    free();
    FILE *nl;
    efunc *r_ops_int[N_OPS];

    //Initialize the AMPL library
    asl = ASL_alloc(ASL_read_fg);

    //Initialize the nl file reading
    nl = jac0dim(const_cast<char*>(stub.c_str()), (fint ) stub.length());

    //Prepare *columnwise* parsing of nl file
    A_vals = (real *) Malloc(nzc * sizeof(real));

    // allocate initial values for primal and dual variables if available
    want_xpi0 = 3;

    // allocate space for initial values
    X0 = new real[n_var];
    havex0 = new char[n_var];
    pi0 = new real[n_con];
    havepi0 = new char[n_con];

#ifdef AMPLDEBUG
    cout << "number of nonzeros    = " << nzc << endl;
    cout << "number of variables   = " << n_var << endl;
    cout << "number of constraints = " << n_con << endl;
    cout << "number of objectives  = " << n_obj << endl;
    cout << "number of ranges      = " << nranges << endl;
    cout << "number of equations   = " << n_eqn << endl;
#endif
    if (N_OPS > 0) {
      for (int i = 0; i < N_OPS; i++) {
        r_ops_int[i] = (efunc*) (unsigned long) i;
      }
      R_OPS = r_ops_int;
      want_derivs = 0;
      fg_read(nl, ASL_keep_all_suffixes);
      R_OPS = 0;
    }
    _isCW = true;
  }
}
void NlReader::getRW() {
  if (_isCW) {
    free();
    FILE *nl;
    efunc *r_ops_int[N_OPS];

    //Initialize the AMPL library
    asl = ASL_alloc(ASL_read_fg);

    //Initialize the nl file reading
    nl = jac0dim(const_cast<char*>(stub.c_str()), (fint ) stub.length());

    //Prepare *columnwise* parsing of nl file
    A_vals = (real *) Malloc(nzc * sizeof(real));

    // allocate initial values for primal and dual variables if available
    want_xpi0 = 3;

    // allocate space for initial values
    X0 = new real[n_var];
    havex0 = new char[n_var];
    pi0 = new real[n_con];
    havepi0 = new char[n_con];

#ifdef AMPLDEBUG
    cout << "number of nonzeros    = " << nzc << endl;
    cout << "number of variables   = " << n_var << endl;
    cout << "number of constraints = " << n_con << endl;
    cout << "number of objectives  = " << n_obj << endl;
    cout << "number of ranges      = " << nranges << endl;
    cout << "number of equations   = " << n_eqn << endl;
#endif
    if (N_OPS > 0) {
      for (int i = 0; i < N_OPS; i++) {
        r_ops_int[i] = (efunc*) (unsigned long) i;
      }
      R_OPS = r_ops_int;
      want_derivs = 0;
      qp_read(nl, ASL_keep_all_suffixes);
      R_OPS = 0;
    }
    _isCW = false;
  }
}

NlReader::NlReader(std::string nlfilename, std::string osol)
    : stub(nlfilename),
      _isCW(false) {
  BuildOppName();
  og = 0;
  asl = 0;
  getCW();
  _onoff = false;

  //	// Now create row-wise version
  //	asl = rw = ASL_alloc(ASL_read_fg);
  //	nl = jac0dim((char*) stub.c_str(), (fint) stub.length());
  //	want_derivs = 0;
  //	qp_read(nl, 0);
  //
  //	asl = cw;
  numkount = 0;
  //	asl = rw;

  // store the osol string into the OSnl2OS object
  this->osol = osol;

}

NlReader::~NlReader() {
  //	if (osinstance != 0)
  //		delete osinstance;
  //	osinstance = 0;
  //
  //	if (osolreader != 0) {
  //		delete osolreader;
  //		osolreader = 0;
  //	} else {
  //		if (osoption != 0)
  //			delete osoption;
  //		osoption = 0;
  //	}
  free();
  //	free(A_vals);
  //	if (X0)
  //		delete[] X0;
  //	X0 = 0;
  //	if (havex0)
  //		delete[] havex0;
  //	havex0 = 0;
  //	if (pi0)
  //		delete[] pi0;
  //	pi0 = 0;
  //	if (havepi0)
  //		delete[] havepi0;
  //	havepi0 = 0;
  //	ASL_free (&cw);
  //	ASL_free(&rw);
}

//1
//2
//16
//46
//54
//76
//79
//81

void NlReader::BuildOppName() {
  if (OppName.empty()) {
    OppName[0] = "OPPLUS";
    OppName[1] = "OPMINUS";
    OppName[2] = "OPMULT";
    OppName[3] = "OPDIV";
    OppName[4] = "OPREM";
    OppName[5] = "OPPOW";
    OppName[6] = "OPLESS";
    OppName[11] = "MINLIST";
    OppName[12] = "MAXLIST";
    OppName[13] = "FLOOR";
    OppName[14] = "CEIL";
    OppName[15] = "ABS";
    OppName[16] = "OPUMINUS";
    OppName[21] = "OPAND";
    OppName[22] = "LT";
    OppName[23] = "LE";
    OppName[24] = "EQ";
    OppName[28] = "GE";
    OppName[29] = "GT";
    OppName[30] = "NE";
    OppName[34] = "OPNOT";
    OppName[35] = "OPIFnl";
    OppName[28] = "GE";
    OppName[29] = "GT";
    OppName[30] = "NE";
    OppName[34] = "OPNOT";
    OppName[35] = "OPIFnl";
    OppName[37] = "OP_tanh";
    OppName[38] = "OP_sqrt";
    OppName[39] = "OP_sinh";
    OppName[40] = "OP_sin";
    OppName[41] = "OP_sin";
    OppName[42] = "OP_log10";
    OppName[43] = "OP_log";
    OppName[44] = "OP_exp";
    OppName[45] = "OP_cosh";
    OppName[46] = "OP_cos";
    OppName[47] = "OP_atanh";
    OppName[48] = "OP_atan2";
    OppName[49] = "OP_atan";
    OppName[50] = "OP_asinh";
    OppName[51] = "OP_asin";
    OppName[52] = "OP_acosh";
    OppName[53] = "OP_acos";
    OppName[54] = "OPSUMLIST";
    OppName[55] = "OPintDIV";
    OppName[56] = "OPprecision";
    OppName[57] = "OPround";
    OppName[58] = "OPtrunc";
    OppName[59] = "OPCOUNT";
    OppName[60] = "OPNUMBEROF";
    OppName[61] = "OPNUMBEROFs";
    OppName[62] = "OPATLEAST";
    OppName[63] = "OPATMOST";
    OppName[64] = "OPPLTERM";
    OppName[65] = "OPIFSYM";
    OppName[66] = "OPEXACTLY";
    OppName[67] = "OPNOTATLEAST";
    OppName[68] = "OPNOTATMOST";
    OppName[69] = "OPNOTEXACTLY";
    OppName[70] = "ANDLIST";
    OppName[71] = "ORLIST";
    OppName[72] = "OPIMPELSE";
    OppName[73] = "OP_IFF";
    OppName[74] = "OPALLDIFF";
    OppName[75] = "OP1POW";
    OppName[76] = "OP2POW";
    OppName[77] = "OPCPOW";
    OppName[78] = "OPFUNCALL";
    OppName[79] = "OPNUM";
    OppName[80] = "OPHOL";
    OppName[81] = "OPVARVAL";
    OppName[82] = "N_OPS";
  }
}

Function NlReader::walkTree(Problem & p, expr *e) {
  efunc *op;
  expr **ep;
  int opnum;
  int j = ((expr_v *) e - VAR_E) - p.variables().size();
  op = e->op;
  opnum = Intcast op;

  Function node;
  //  std::cout << "in  node : " << node << std::endl;
  //  std::cout << "operator : "<<OppName[OPPLUS]<<std::endl;
  switch (opnum) {
    case OPPLUS:
      node = walkTree(p, e->L.e);
      node += walkTree(p, e->R.e);
      break;
    case OPSUMLIST:
      for (ep = e->L.ep; ep < e->R.ep; ep++) {
        //			for (ep = e->L.ep; ep < e->R.ep; *ep++) {
        node += walkTree(p, *ep);
      }
      break;
    case MAXLIST:
      UN_IMPLEMENTTED_OPP(MAXLIST);
      break;
    case OPMINUS:
      node = walkTree(p, e->L.e);
      node -= walkTree(p, e->R.e);
      break;
    case OPUMINUS:
      node = walkTree(p, e->L.e);
      node = -node;
      break;
    case OPMULT:
      node = walkTree(p, e->L.e);
      node *= walkTree(p, e->R.e);
      break;
    case OPDIV:
      UN_IMPLEMENTTED_OPP(OPDIV);
      //      node = walkTree(p, e->L.e);
      //      node /= walkTree(p, e->R.e);
      break;
    case OPPOW:
      UN_IMPLEMENTTED_OPP(OPPPOW);
      break;
    case OP1POW:
      UN_IMPLEMENTTED_OPP(OP1POW);
      //      node = walkTree(p, e->L.e);
      //      node = pow(node, e->R.en->v);
      break;
    case OP2POW:
      UN_IMPLEMENTTED_OPP(OP2POW);
      //      node = walkTree(p, e->L.e);
      //      node *= walkTree(p, e->L.e);
      break;
    case OP_log:
      UN_IMPLEMENTTED_OPP(OP_log);
      //      node = walkTree(p, e->L.e);
      //      node = log(node);
      break;
    case OP_sqrt:
      UN_IMPLEMENTTED_OPP(OP_sqrt);
      //      node = walkTree(p, e->L.e);
      //      node = pow(node, 0.5);
      //		std::cerr << "NOT IMPLEMENTED OP_sqrt" << std::endl;
      break;
    case OP_cos:
      UN_IMPLEMENTTED_OPP(OP_cos);
      //      node = walkTree(p, e->L.e);
      //      node = cos(node);
      break;
    case OP_sin:
      UN_IMPLEMENTTED_OPP(OP_sin);
      //      node = walkTree(p, e->L.e);
      //      node = sin(node);
      break;
    case OP_exp:
      UN_IMPLEMENTTED_OPP(OP_exp);
      //      node = walkTree(p, e->L.e);
      //      node = exp(node);
      break;
    case ABS:
      UN_IMPLEMENTTED_OPP(ABS);
      break;
    case OPNUM:
      node += (double) ((expr_n*) e)->v;
      break;
    case OPVARVAL:
      if (j >= 0) {
        // process common expression

        //						std::cout << "como = " << como << std::endl;
        //						std::cout << "comc = " << comc << std::endl;
        //						std::cout << "comb = " << comb << std::endl;
        //						std::cout << "como1 = " << como1 << std::endl;
        //						std::cout << "comc1 = " << comc1 << std::endl;
        //						std::cout << "ncom0 = " << ncom0 << std::endl;
        //						std::cout << "jjjjjjjjjjjjjjjjjj = " << j << std::endl;

        // Orban: http://www.gerad.ca/~orban/drampl/def-vars.html
        if (j < ncom0) {
          struct cexp *common = ((const ASL_fg *) asl)->I.cexps_ + j;
          int nlin = common->nlin;

          if (nlin > 0) {
            linpart *L = common->L;
            for (int kj = 0; kj < nlin; kj++) {
              size_t const idx(
                  ((uintptr_t)(L[kj].v.rp) - (uintptr_t) VAR_E)
                      / sizeof(expr_v));
              double const coeff(L[kj].fac);
              node += coeff * p.variable(idx);
            }
            node += walkTree(p, common->e);
          } else {
            node += walkTree(p, common->e);
          }
        } else {
          struct cexp1 *common = ((const ASL_fg *) asl)->I.cexps1_
              + (j - ncom0);
          int nlin = common->nlin;
          if (nlin > 0) {
            linpart *L = common->L;
            for (int kj = 0; kj < nlin; kj++) {
              size_t const idx(
                  ((uintptr_t)(L[kj].v.rp) - (uintptr_t) VAR_E)
                      / sizeof(expr_v));
              double
              const coeff(L[kj].fac);
              node += coeff * p.variable(idx);
            }
            node += walkTree(p, common->e);
          } else {
            node += walkTree(p, common->e);
          }
        }
      } else {
        assert(e->a < (int ) p.variables().size());
        node += p.variable(e->a);
      }
      break;
    default:
      std::cout
          << "ERROR:  An unsupported operator found, AMPL operator number =  "
          << opnum << std::endl;
      break;  //end switch
  }
  //  std::cout << "out node : " << node << std::endl;
  //	if (_onoff)
  //		std::cout << "OPNUM " << OppName[opnum] << " : " << node << std::endl;
  return node;
}                            //walkTree

//static inline char integerOrBinary(real upper, real lower) {
//	if (lower > -1.0 + MyEps && upper < 2.0 - MyEps)
//		return 'B';
//	return 'I';
//}

bool NlReader::createVariables(Problem & p) {
  p.variables().reserve(n_var);
  //    osinstance->setVariableNumber( n_var);

  //first the nonlinear variables
  //welcome to the world of the ASL API
  int i;
  int lower;
  int upper;
  lower = 0;
  upper = nlvb - nlvbi;
#ifdef AMPLDEBUG
  std::cout << "LOWER = " << lower << std::endl;
  std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //continuous in an objective and in a constraint
      {
    //		std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
    //        osinstance->addVariable(i, var_name(i),
    //                                LUv[2*i]   > -MyInfinity ? LUv[2*i]   : -MyInfinity,
    //                                LUv[2*i+1] <  MyInfinity ? LUv[2*i+1] :  MyInfinity,
    //                                // vartype);
  }

  lower = nlvb - nlvbi;
  upper = (nlvb - nlvbi) + nlvbi;
  upper = nlvb;  //
#ifdef AMPLDEBUG
      std::cout << "LOWER = " << lower << std::endl;
      std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //integer in an objective and in a constraint
      {

    // vartype = integerOrBinary(LUv[2 * i], LUv[2 * i + 1]); // AMPL doesn't make the distinction for nonlinear variables
    //		osinstance->addVariable(i, var_name(i),
    //				LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    //				LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity,
    //				// vartype);
    //		std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
  }

  lower = nlvb;
  upper = nlvb + (nlvc - (nlvb + nlvci));
  upper = nlvc - nlvci;
#ifdef AMPLDEBUG
  std::cout << "LOWER = " << lower << std::endl;
  std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //continuous just in constraints
      {

    // vartype = 'C';
    //		osinstance->addVariable(i, var_name(i),
    //				LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    //				LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity,
    //				// vartype);
    std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
  }

  lower = nlvc - nlvci;
  upper = nlvc - nlvci + nlvci;
  upper = nlvc;
#ifdef AMPLDEBUG
  std::cout << "LOWER = " << lower << std::endl;
  std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //integer just in constraints
      {

    // vartype = integerOrBinary(LUv[2 * i], LUv[2 * i + 1]); // AMPL doesn't make the distinction for nonlinear variables
    //		osinstance->addVariable(i, var_name(i),
    //				LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    //				LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity,
    //				// vartype);
    std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
  }

  lower = nlvc;
  upper = nlvc + (nlvo - (nlvc + nlvoi));
  upper = nlvo - nlvoi;
#ifdef AMPLDEBUG
  std::cout << "LOWER = " << lower << std::endl;
  std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //continuous just in objectives
      {

    // vartype = 'C';
    //		osinstance->addVariable(i, var_name(i),
    //				LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    //				LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity,
    //				// vartype);
    std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
  }

  lower = nlvo - nlvoi;
  upper = nlvo - nlvoi + nlvoi;
  upper = nlvo;
#ifdef AMPLDEBUG
  std::cout << "LOWER = " << lower << std::endl;
  std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //integer just in objectives
      {

    // vartype = integerOrBinary(LUv[2 * i], LUv[2 * i + 1]); // AMPL doesn't make the distinction for nonlinear variables
    //		osinstance->addVariable(i, var_name(i),
    //				LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    //				LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity,
    //				// vartype);

    std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
  }

  //now the other variables

  lower = std::max(nlvc, nlvo);
  upper = std::max(nlvc, nlvo) + nwv;
#ifdef AMPLDEBUG
  std::cout << "LOWER = " << lower << std::endl;
  std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //linear arc variables
      {

    // vartype = 'C';
    //		osinstance->addVariable(i, var_name(i),
    //				LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    //				LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity,
    //				// vartype);
    std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
  }

  lower = std::max(nlvc, nlvo) + nwv;
  upper = std::max(nlvc, nlvo) + nwv
      + (n_var - (std::max(nlvc, nlvo) + niv + nbv + nwv));
  upper = n_var - niv - nbv;
#ifdef AMPLDEBUG
  std::cout << "LOWER = " << lower << std::endl;
  std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //other linear
      {

    // vartype = 'C';
    //		osinstance->addVariable(i, var_name(i),
    //				LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    //				LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity,
    //				// vartype);
    std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
  }

  lower = n_var - niv - nbv;
  upper = n_var - niv - nbv + nbv;
  upper = n_var - niv;
#ifdef AMPLDEBUG
  std::cout << "LOWER = " << lower << std::endl;
  std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //linear binary
      {

    // vartype = 'B';
    //		osinstance->addVariable(i, var_name(i),
    //				LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    //				LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity,
    //				// vartype);
    std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
  }

  lower = n_var - niv;
  upper = n_var - niv + niv;
  upper = n_var;
#ifdef AMPLDEBUG
  std::cout << "LOWER = " << lower << std::endl;
  std::cout << "UPPER = " << upper << std::endl;
#endif
  for (i = lower; i < upper; i++)  //linear integer
      {

    // vartype = 'I';
    //		osinstance->addVariable(i, var_name(i),
    //				LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    //				LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity,
    //				// vartype);

    //		std::string s(var_name(i));
    p.newVariable(LUv[2 * i] > -MyInfinity ? LUv[2 * i] : -MyInfinity,
    LUv[2 * i + 1] < MyInfinity ? LUv[2 * i + 1] : MyInfinity);
  }

  // end of variables -- thank goodness!!!
  return true;
}

bool NlReader::createObjectif(Problem & p) {
  getRW();
  // now create the objective function
  // in the nl file, this is stored in dense form; convert to sparse.
  //
  std::vector<Function> builder;
  builder.reserve(n_obj);
//  std::cout << "n_obj = " << n_obj << std::endl;
  for (int i(0); i < n_obj && i == 0; i++) {
    builder.push_back(walkTree(p, OBJ_DE[i].e));
    for (og = Ograd[i]; og; og = og->next) {
      if (fabs(og->coef) > MyEps) {
        double const v((objtype[i] == 1) ? -og->coef : og->coef);
        builder[i] += p.variable(og->varno) * v;
      }
    }
  }

  p.objectifs().clear();
  for (int i(0); i < n_obj && i == 0; i++) {
    //		std::cout << "o" << i << " : " << builder[i] << std::endl;
    assert(i == (int ) p.addObj(builder[i]));
  }
  return true;
}
bool NlReader::createConstraints(Problem &p) {
  //Switch to row-wise format.
  //	asl = rw;
  getRW();
  //
  // now fill in row information
  //
  std::vector<Function> builder;
  builder.reserve(n_con);
  //	p.constraints().reserve(n_con);
  //		osinstance->setConstraintNumber(n_con);
  // kipp -- important  -- figure out where the nl file stores a rhs constant

  for (int i(0); i < n_con; i++) {
    builder.push_back(walkTree(p, CON_DE[i].e));
//    std::cout << "builder[" << i << "] : " << builder[i] << std::endl;
    //		if (i == 9)
    //			exit(0);
    //    		std::cout << builder[i] << std::endl;
  }

  //	j = A_colstarts[i-1] - 1;
  //	for(je = A_colstarts[i] - 1; j < je; j++)
  //		printf("\t%ld\t%g\n", A_rownos[j], A_vals[j]);
  getCW();
  for (int i(0); i < n_var; i++) {
    for (int j(A_colstarts[i]); j < A_colstarts[i + 1]; j++) {
      builder[A_rownos[j]] += (A_vals[j] * p.variable(i));
    }
  }
  getRW();
  for (int i(0); i < n_con; i++) {
    //    		std::cout << builder[i] << std::endl;
    //		std::cout << "CONSTRAINT " << i << std::endl;
    //		std::string const name(con_name(i));
    Constraint constraint(
        builder[i], LUrhs[2 * i] > -MyInfinity ? LUrhs[2 * i] : -MyInfinity,
        LUrhs[2 * i + 1] < MyInfinity ? LUrhs[2 * i + 1] : MyInfinity);
    p.newConstraint(constraint);
    //		std::cout << constraint << std::endl;
    //		if (i == 9)
    //			exit(0);
  }
  //	exit(0);
  //	std::cout << "CONSTRAINT " << 22534 << " : " << std::endl;
  //	std::cout << builder[22534] << std::endl;
  return true;
}
bool NlReader::createObjects(Problem & p) {
  p.clear();
  createVariables(p);
  createObjectif(p);
  createConstraints(p);
  // process initial primal values

  // count the number of values (including those in the OSoL file)
  int n_x0(0);
  for (int i(0); i < n_var; i++)
    if (havex0[i] != 0)
      n_x0++;
  //  p.initPrimal().clear();
  //  if (n_x0 > 0) {
  //    // pull info out of ASL data structure
  //    for (int i(0); i < n_var; i++)
  //      if (havex0[i] != 0) {
  //        p.initPrimal()[i] = X0[i];
  //      }
  //  }

  // count the number of values (including those in the OSoL file)
  int n_pi0 = 0;
  for (int i = 0; i < n_con; i++)
    if (havepi0[i] != 0)
      n_pi0++;

  //  if (n_pi0 > 0) {
  //    for (int i = 0; i < n_con; i++) {
  //      if (havepi0[i] != 0) {
  //        p.initDual()[i] = pi0[i];
  //      }
  //    }
  //  }
  //
  //	/** Before we can process the rest of the instance, we check for QP
  //	 *  This needs to be done here because of the possibility of expressions like (1 - x[0])^2
  //	 *  which may modify the A matrix as well as the right-hand sides.
  //	 *  If the A-matrix is modified, the column-wise representation
  //	 *  is out of date and must be rebuilt from the row-wise form.
  //	 */
  //	std::vector<int> fidxs, v1idxs, v2idxs;
  //	std::vector<double> coeffs;
  ////		std::vector<Nl> nlExprs;
  //	std::vector<Tree2> nlExprs;
  //	real* delsqp;
  //	fint* colqp;
  //	fint* rowqp;
  //	int osNLIdx; // OS n.l. function index
  //	int aNLIdx; // AMPL n.l. function index
  //
  ////Switch to row-wise format.
  ////	asl = rw;
  //	getRW();
  //
  //// Iterate from -nlo to nlc-1 so that the qterms are sorted by idx
  ////	{
  //	expr* e = CON_DE [24].e; // because osNLIdx = -aNLIdx-1
  //
  ////	_onoff = true;
  ////	walkTree(p, e);
  ////	exit(0);
  ////	_onoff = false;
  ////	}
  //// Process the objectives first, for which fill-in does not matter
  //	for (osNLIdx = -nlo, aNLIdx = nlo - 1; osNLIdx < 0; osNLIdx++, aNLIdx--) {
  //		//		if (nqpcheck(aNLIdx, &rowqp, &colqp, &delsqp) > 0) // quadratic
  //		//				{
  //		//			for (int v1 = 0; v1 < n_var; v1++) {
  //		//				for (int* psV2 = &rowqp[colqp[v1]];
  //		//						psV2 < &rowqp[colqp[v1 + 1]]; psV2++, delsqp++) {
  //		//					if (std::abs(*delsqp) > MyEps) // Try to exclude terms introduced by rounding
  //		//							{
  //		//						fidxs.push_back(osNLIdx);
  //		//						v1idxs.push_back(v1);
  //		//						v2idxs.push_back(*psV2);
  //		//						coeffs.push_back(0.5 * *delsqp);
  //		//					}
  //		//				}
  //		//			}
  //		//		} else // Nonlinear or error in nqpcheck
  //		//		{
  //		expr* e = aNLIdx < 0 ? CON_DE [osNLIdx].e : OBJ_DE [aNLIdx].e; // because osNLIdx = -aNLIdx-1
  //		//			nl.idx = osNLIdx;
  //		//			nl.osExpressionTree = new OSExpressionTree();
  //		//			nl.osExpressionTree->m_treeRoot = walkTree(e);
  //		//			nl.m_bDeleteExpressionTree = false;
  //
  //		/*
  //		 * Note: If the copy operation of the Nl class is changed from shallow
  //		 * to deep, we will want to manage memory differently here.
  //		 */
  //		nlExprs.push_back(walkTree(p, e));
  //		//		}
  //	}
  //
  //	bool isQP = true;
  //	bool fill_in = false;
  //	int nqpchk;
  //	cgrad *cg;
  //
  //	double* A_row_temp = new double[n_var];
  //
  //	for (osNLIdx = 0, aNLIdx = -1; osNLIdx < nlc; osNLIdx++, aNLIdx--) {
  //		//		if (isQP) // No need to identify quadratic terms once we have found a non-quadratic term
  //		//		{
  //		//			// check the nonzeroes before and after
  //		//			if (!fill_in) // once we know there will be fill-in we can stop counting
  //		//			{
  //		//				for (cg = Cgrad[osNLIdx]; cg; cg = cg->next) {
  //		//					if (cg->coef != 0)
  //		//						A_row_temp[cg->varno] = cg->coef;
  //		//				}
  //		//			}
  //		//
  //		//			nqpchk = nqpcheck(aNLIdx, &rowqp, &colqp, &delsqp);
  //		//			if (nqpchk > 0) // quadratic
  //		//					{
  //		//				for (int v1 = 0; v1 < n_var; v1++) {
  //		//					for (int* psV2 = &rowqp[colqp[v1]];
  //		//							psV2 fill_in< &rowqp[colqp[v1 + 1]]; psV2++, delsqp++) {
  //		//						if (std::abs(*delsqp) > MyEps) // Try to exclude terms introduced by rounding
  //		//								{
  //		//							fidxs.push_back(osNLIdx);
  //		//							v1idxs.push_back(v1);
  //		//							v2idxs.push_back(*psV2);
  //		//							coeffs.push_back(0.5 * *delsqp);
  //		//						}
  //		//					}
  //		//				}
  //		//				if (!fill_in) // once we know there will be fill-in we can stop counting
  //		//				{
  //		//					for (cg = Cgrad[osNLIdx]; cg; cg = cg->next) {
  //		//						if (cg->coef != 0)
  //		//							if (cg->coef != A_row_temp[cg->varno]) {
  //		//								fill_in = true;
  //		//								break;
  //		//							}
  //		//					}
  //		//				}
  //		//				continue;
  //		//			}
  //		//			if (nqpchk < 0)
  //		//				isQP = false;
  //		//		}
  //
  //		// Nonlinear or error in nqpcheck
  //		{
  //			//			Nl nl;
  //			expr* e = aNLIdx < 0 ? CON_DE [osNLIdx].e : OBJ_DE [aNLIdx].e; // because osNLIdx = -aNLIdx-1
  //			//			nl.idx = osNLIdx;
  //			//			nl.osExpressionTree = new OSExpressionTree();
  //			//			nl.osExpressionTree->m_treeRoot = walkTree(e);
  //			//			nl.m_bDeleteExpressionTree = false;
  //			//			/*
  //			//			 * Note: If the copy operation of the Nl class is changed from shallow
  //			//			 * to deep, we will want to manage memory differently here.
  //			//			 */
  //			//			nlExprs.push_back(nl);
  //			nlExprs.push_back(walkTree(p, e));
  //		}
  //	}
  //	delete[] A_row_temp;
  //
  ////	if (nlExprs.size()) {
  ////		Nl** ppsNl = new Nl*[nlExprs.size()];
  ////		for (i = 0; i < nlExprs.size(); i++) {
  ////			ppsNl[i] = new Nl(nlExprs[i]); // See above note about shallow copy
  ////			ppsNl[i]->m_bDeleteExpressionTree = true;
  ////		}
  ////		osinstance->instanceData->nonlinearExpressions->nl = ppsNl;
  ////	}
  ////	//	osinstance->instanceData->nonlinearExpressions->numberOfNonlinearExpressions =
  ////	//			nlExprs.size();
  ////	//	if (fidxs.size()) {
  ////	//		osinstance->setQuadraticTerms((int) fidxs.size(), &fidxs[0], &v1idxs[0],
  ////	//				&v2idxs[0], &coeffs[0], 0, (int) fidxs.size() - 1);
  ////	//	}
  ////	//	// Note: if we intended to call objval, conval etc with asl == rw later we must call qp_opify here.
  ////	//
  ////	//	//
  ////	//	// end loop of nonlinear rows
  ////	//	//
  //
  //// now create the objective function
  //// in the nl file, this is stored in dense form; convert to sparse.
  ////
  //
  ////	for (i = 0; i < n_obj && i == 0; i++) {
  //	IdxType const nObj(0);
  //	p._minimize = walkTree(p, OBJ_DE [nObj].e);
  //	for (og = Ograd[nObj]; og; og = og->next) {
  //		if (fabs(og->coef) > MyEps) {
  //			p._minimize += p.variable(og->varno)
  //					* (og->coef * (objtype[i] == 1) ? -1 : 1);
  //		}
  //	}
  ////	}
  //
  ////
  //// now fill in row information
  ////
  //	p.constraints().reserve(n_con);
  ////		osinstance->setConstraintNumber(n_con);
  //// kipp -- important  -- figure out where the nl file stores a rhs constant
  //	for (i = 0; i < n_con; i++) {
  //		std::string const name(con_name(i));
  //		Tree2 xp(walkTree(p, CON_DE [i].e));
  ////		for (cg = Cgrad[i]; cg; cg = cg->next) {
  //////			if (cg->coef != 0) {
  //////				xp += p.variable(cg->varno) * cg->coef;
  //////			}
  ////		}
  //		Constraint constraint(xp,
  //				LUrhs[2 * i] > -MyInfinity ? LUrhs[2 * i] : -MyInfinity,
  //				LUrhs[2 * i + 1] < MyInfinity ? LUrhs[2 * i + 1] : MyInfinity);
  //		p.newConstraint(constraint);
  //	}
  //
  //	int row_len;
  ////	A_rowstarts = new int[n_con + 1];
  ////	A_rowstarts[0] = 0;
  ////	for (int i = 0; i < n_con; i++) {
  ////		row_len = 0;
  ////		for (cg = Cgrad[i]; cg; cg = cg->next) {
  ////			if (cg->coef != 0)
  ////				row_len++;
  ////		}
  ////		A_rowstarts[i + 1] = A_rowstarts[i] + row_len;
  ////	}
  ////	for (int i = 0; i < n_con; i++) {
  ////		for (cg = Cgrad[i]; cg; cg = cg->next) {
  ////			if (cg->coef != 0) {
  ////				p.constraint(i) += cg->coef * p.variable(cg->varno);
  ////			}
  ////		}
  ////	}
  //	getCW();
  //	int colEnd = 0;
  //	for (i = 0; i < n_var; i++) {
  //		int colStart = colEnd;
  //		int colEnd = A_colstarts[i + 1];
  //		for (j = colStart; j < colEnd; j++) {
  //			if (fabs(A_vals[j]) > MyEps) {
  //				p.constraint(A_rownos[j]) += A_vals[j] * p.variable(i);
  ////				exit(0);
  //			}
  //		}
  //	}
  ////	//
  ////	// Now the A-matrix
  ////	// The treatment depends on whether there was fill-in during the QP check or not
  ////	//
  ////	if (fill_in) // store the matrix rowwise
  ////	{
  ////		int row_len;
  ////		A_rowstarts = new int[n_con + 1];
  ////		A_rowstarts[0] = 0;
  ////		for (int i = 0; i < n_con; i++) {
  ////			row_len = 0;
  ////			for (cg = Cgrad[i]; cg; cg = cg->next) {
  ////				if (cg->coef != 0)
  ////					row_len++;
  ////			}
  ////			A_rowstarts[i + 1] = A_rowstarts[i] + row_len;
  ////		}
  ////		A_colptr = new int[A_rowstarts[n_con]];
  ////		A_nzelem = new double[A_rowstarts[n_con]];
  ////		for (int i = 0; i < n_con; i++) {
  ////			row_len = 0;
  ////			for (cg = Cgrad[i]; cg; cg = cg->next) {
  ////				if (cg->coef != 0) {
  ////					A_colptr[A_rowstarts[i] + row_len] = cg->varno;
  ////					A_nzelem[A_rowstarts[i] + row_len] = cg->coef;
  ////					row_len++;
  ////				}
  ////			}
  ////		}
  ////
  ////		if (A_rowstarts[n_con] > 0) {
  ////			osinstance->setLinearConstraintCoefficients(A_rowstarts[n_con],
  ////					false, A_nzelem, 0, A_rowstarts[n_con] - 1, A_colptr, 0,
  ////					A_rowstarts[n_con] - 1, A_rowstarts, 0, n_con);
  ////			// setLinearConstraintCoefficients does a soft copy, and unlike the column-wise representation,
  ////			// this row-wise representation is not taken care of by the ASL_free method in the destructor
  ////			osinstance->instanceData->linearConstraintCoefficients->start->bDeleteArrays =
  ////					true;
  ////			osinstance->instanceData->linearConstraintCoefficients->colIdx->bDeleteArrays =
  ////					true;
  ////			osinstance->instanceData->linearConstraintCoefficients->value->bDeleteArrays =
  ////					true;
  ////		}
  ////
  ////#ifdef AMPLDEBUG
  ////		cout << "A-matrix elements: ";
  ////		for (int i = 0; i < A_rowstarts[n_con]; i++)
  ////			cout << A_nzelem[i] << " ";
  ////		cout << endl;
  ////		cout << "A-matrix col index: ";
  ////		for (int i = 0; i < A_rowstarts[n_con]; i++)
  ////			cout << A_colptr[i] << " ";
  ////		cout << endl;
  ////		cout << "A-matrix rowstart: ";
  ////		for (int i = 0; i <= n_con; i++)
  ////			cout << A_rowstarts[i] << " ";
  ////		cout << endl;
  ////#endif
  ////}
  ////
  ////	else {
  ////		asl = cw;
  ////		int colStart, colEnd, nCoefSqueezed;
  ////		nCoefSqueezed = 0;
  ////
  ////	#ifdef AMPLDEBUG
  ////			cout << "A-matrix elements: ";
  ////			for (int i = 0; i < A_colstarts[n_var]; i++)
  ////				cout << A_vals[i] << " ";
  ////			cout << endl;
  ////			cout << "A-matrix rowinfo: ";
  ////			for (int i = 0; i < A_colstarts[n_var]; i++)
  ////				cout << A_rownos[i] << " ";
  ////			cout << endl;
  ////			cout << "A-matrix colstart: ";
  ////			for (int i = 0; i <= n_var; i++)
  ////				cout << A_colstarts[i] << " ";
  ////			cout << endl;
  ////	#endif
  ////
  ////		colEnd = 0;
  ////		for (i = 0; i < n_var; i++) {
  ////			colStart = colEnd;
  ////			colEnd = A_colstarts[i + 1];
  ////#ifdef AMPLDEBUG
  ////			cout << "col " << i << " from " << colStart << " to " << colEnd - 1
  ////					<< endl;
  ////#endif
  ////			for (j = colStart; j < colEnd; j++) {
  ////				if (fabs(A_vals[j]) > MyEps) {
  ////					A_vals[j - nCoefSqueezed] = A_vals[j];
  ////					A_rownos[j - nCoefSqueezed] = A_rownos[j];
  ////				} else {
  ////#ifdef AMPLDEBUG
  ////					cout << "squeeze out element " << j << endl;
  ////#endif
  ////					nCoefSqueezed++;
  ////				}
  ////			}
  ////			A_colstarts[i + 1] = A_colstarts[i + 1] - nCoefSqueezed;
  ////		}
  ////
  ////#ifdef AMPLDEBUG
  ////		cout << "A-matrix elements: ";
  ////		for (i = 0; i < A_colstarts[n_var]; i++)
  ////			cout << A_vals[i] << " ";
  ////		cout << endl;
  ////		cout << "A-matrix rowinfo: ";
  ////		for (i = 0; i < A_colstarts[n_var]; i++)
  ////			cout << A_rownos[i] << " ";
  ////		cout << endl;
  ////		cout << "A-matrix colstart: ";
  ////		for (i = 0; i <= n_var; i++)
  ////			cout << A_colstarts[i] << " ";
  ////		cout << endl;
  ////		cout << "A-matrix nonzeroes: " << A_colstarts[n_var] << "; nsqueezed: "
  ////				<< nCoefSqueezed << endl;
  ////#endif
  ////
  ////		if (A_colstarts[n_var] > 0) {
  ////			osinstance->setLinearConstraintCoefficients(A_colstarts[n_var],
  ////					true, A_vals, 0, A_colstarts[n_var] - 1, A_rownos, 0,
  ////					A_colstarts[n_var] - 1, A_colstarts, 0, n_var);
  ////			osinstance->instanceData->linearConstraintCoefficients->start->bDeleteArrays =
  ////					false;
  ////			osinstance->instanceData->linearConstraintCoefficients->rowIdx->bDeleteArrays =
  ////					false;
  ////			osinstance->instanceData->linearConstraintCoefficients->value->bDeleteArrays =
  ////					false;
  ////		}
  ////	}
  ////
  ////	/**
  ////	 *  The nl file may contain options that are indexed over model entities:
  ////	 *  variables, constraints, objectives, problems. An example would be
  ////	 *  initial primal and dual variables. AMPL stores these as suffixes,
  ////	 *  which are processed next. The code below is based on ideas expressed by David Gay.
  ////	 */
  ////
  //	SufDesc *d;
  //	int suffixType, nOther, nOtherIdx;
  //
  ////		asl = cw;
  //	getCW();
  //
  ////	try {
  ////		osolreader = new OSoLReader();
  ////
  //// check if there are any suffixes to deal with and read options if necessary
  //	if ((asl->i.suffixes[ASL_Sufkind_var] != 0)
  //			|| (asl->i.suffixes[ASL_Sufkind_con] != 0)
  //			|| (asl->i.suffixes[ASL_Sufkind_obj] != 0)
  //			|| (asl->i.suffixes[ASL_Sufkind_prob] != 0)) {
  ////				if (osol != "")
  ////					osoption = osolreader->readOSoL(osol);
  ////				else
  ////					osoption = new OSOption();
  //		std::cout << "ya des suffixes" << __FILE__ << " l. " << __LINE__
  //				<< std::endl;
  //	} else {
  //		std::cout << "ya pas de suffixe " << __FILE__ << " l. " << __LINE__
  //				<< std::endl;
  //	}
  ////
  ////		bool found;
  ////		bool extend;
  ////		int nOther;
  ////		std::string *otherOptionNames = 0;
  ////
  ////		// make a record of all <otherVariableOptions> present in the OSoL file
  ////		if (osoption != 0 && osoption->optimization != 0
  ////				&& osoption->optimization->variables != 0
  ////				&& osoption->optimization->variables->numberOfOtherVariableOptions
  ////						> 0) {
  ////			otherOptionNames =
  ////					new std::string[osoption->optimization->variables->numberOfOtherVariableOptions];
  ////			nOther = 0;
  ////			for (int i = 0;
  ////					i
  ////							< osoption->optimization->variables->numberOfOtherVariableOptions;
  ////					i++)
  ////				if (osoption->optimization->variables->other[i]->numberOfVar
  ////						> 0)
  ////					otherOptionNames[nOther++] =
  ////							osoption->optimization->variables->other[i]->name;
  ////		}
  ////
  ////		// First the variable-indexed suffixes
  //	suffixType = ASL_Sufkind_var;
  //	if ((asl->i.suffixes[suffixType] != 0)) {
  //		std::cout << "ya des suffixes" << __FILE__ << " l. " << __LINE__
  //				<< std::endl;
  //		//			OtherVariableOption* varopt;
  //		//			for (d = asl->i.suffixes[suffixType]; d; d = d->next) {
  //		//#ifdef AMPLDEBUG
  //		//				std::cout << "Detected suffix " << d->sufname << "; kind = "
  //		//						<< d->kind << std::endl;
  //		//#endif
  //		//
  //		//				// allocate space
  //		//				varopt = new OtherVariableOption();
  //		//
  //		//				varopt->name = d->sufname;
  //		//				varopt->numberOfEnumerations = 0;
  //		//				varopt->var = new OtherVarOption*[n_var];
  //		//
  //		//				// check if the option was present in the OSoL file
  //		//				found = false;
  //		//				int iopt;
  //		//				for (iopt = 0; iopt < nOther; iopt++) {
  //		//					if (d->sufname == otherOptionNames[iopt]) {
  //		//						found = true;
  //		//						break;
  //		//					}
  //		//				}
  //		//
  //		//				// merge values by overwriting .nl file info
  //		//				if (found) {
  //		//					OtherVariableOption* otherOption;
  //		//					otherOption = osoption->getOtherVariableOption(iopt);
  //		//					for (int i = 0; i < otherOption->numberOfVar; i++) {
  //		//						if (d->kind & 4) // bit-wise mask to distinguish real from integer
  //		//								{
  //		//							d->u.r[otherOption->var[i]->idx] = os_strtod(
  //		//									otherOption->var[i]->value.c_str(), 0);
  //		//						} else
  //		//							d->u.i[otherOption->var[i]->idx] = (int) os_strtod(
  //		//									otherOption->var[i]->value.c_str(), 0)
  //		//									+ 0.1;
  //		//					}
  //		//					if (otherOption->description == "")
  //		//						varopt->description = "combined from osol and .nl data";
  //		//					else
  //		//						varopt->description = otherOption->description
  //		//								+ "; merged with .nl data";
  //		//				} else
  //		//					varopt->description = "transferred from .nl file";
  //		//
  //		//				// count the number of entries
  //		//				if (d->kind & 4) // bit-wise mask to distinguish real from integer
  //		//						{
  //		//					varopt->type = "real";
  //		//					nOtherIdx = 0;
  //		//					for (int k = 0; k < n_var; k++) {
  //		//						if (d->u.r[k] != 0) {
  //		//							varopt->var[nOtherIdx] = new OtherVarOption();
  //		//							varopt->var[nOtherIdx]->idx = k;
  //		//							varopt->var[nOtherIdx]->value = os_dtoa_format(
  //		//									d->u.r[k]);
  //		//							nOtherIdx++;
  //		//						}
  //		//					}
  //		//				} else             // here the suffix values are integer
  //		//				{
  //		//					varopt->type = "integer";
  //		//					nOtherIdx = 0;
  //		//					for (int k = 0; k < n_var; k++) {
  //		//						if (d->u.i[k] != 0) {
  //		//							varopt->var[nOtherIdx] = new OtherVarOption();
  //		//							varopt->var[nOtherIdx]->idx = k;
  //		//							varopt->var[nOtherIdx]->value = os_dtoa_format(
  //		//									(double) d->u.i[k]);
  //		//							nOtherIdx++;
  //		//						}
  //		//					}
  //		//				}
  //		//
  //		//				varopt->numberOfVar = nOtherIdx;
  //		//
  //		//				if (found) {
  //		//					// here we just replace the <var> element and update the numberOfVar
  //		//					delete[] osoption->optimization->variables->other[iopt]->var;
  //		//					osoption->optimization->variables->other[iopt]->var =
  //		//							varopt->var;
  //		//					osoption->optimization->variables->other[iopt]->numberOfVar =
  //		//							nOtherIdx;
  //		//					varopt->var = 0;
  //		//					osoption->optimization->variables->other[iopt]->description =
  //		//							varopt->description;
  //		//				} else {
  //		//					if (!osoption->setAnOtherVariableOption(varopt))
  //		//						throw ErrorClass(
  //		//								"OSnl2OS: Error transfering suffixes on variables");
  //		//				}
  //		//
  //		//				delete varopt;
  //		//				varopt = 0;
  //		//			}
  //	}
  ////
  ////		// suffixes indexed over constraints and objectives work the same way
  //	suffixType = ASL_Sufkind_con;
  //	if ((asl->i.suffixes[suffixType] != 0)) {
  //		std::cout << "ya des suffixes" << __FILE__ << " l. " << __LINE__
  //				<< std::endl;
  //		//			OtherConstraintOption* conopt;
  //		//			for (d = asl->i.suffixes[suffixType]; d; d = d->next) {
  //		//#ifdef AMPLDEBUG
  //		//				std::cout << "Detected suffix " << d->sufname << "; kind = "
  //		//						<< d->kind << std::endl;
  //		//#endif
  //		//
  //		//				// allocate space
  //		//				conopt = new OtherConstraintOption();
  //		//
  //		//				conopt->name = d->sufname;
  //		//				conopt->numberOfEnumerations = 0;
  //		//				conopt->con = new OtherConOption*[n_con];
  //		//
  //		//				// check if the option was present in the OSoL file
  //		//				found = false;
  //		//				int iopt;
  //		//				for (iopt = 0; iopt < nOther; iopt++) {
  //		//					if (d->sufname == otherOptionNames[iopt]) {
  //		//						found = true;
  //		//						break;
  //		//					}
  //		//				}
  //		//
  //		//				// merge values by overwriting .nl file info
  //		//				if (found) {
  //		//					OtherConstraintOption* otherOption;
  //		//					otherOption = osoption->getOtherConstraintOption(iopt);
  //		//					for (int i = 0; i < otherOption->numberOfCon; i++) {
  //		//						if (d->kind & 4) // bit-wise mask to distinguish real from integer
  //		//								{
  //		//							d->u.r[otherOption->con[i]->idx] = os_strtod(
  //		//									otherOption->con[i]->value.c_str(), 0);
  //		//						} else
  //		//							d->u.i[otherOption->con[i]->idx] = (int) os_strtod(
  //		//									otherOption->con[i]->value.c_str(), 0)
  //		//									+ 0.1;
  //		//					}
  //		//					if (otherOption->description == "")
  //		//						conopt->description = "combined from osol and .nl data";
  //		//					else
  //		//						conopt->description = otherOption->description
  //		//								+ "; merged with .nl data";
  //		//				} else
  //		//					conopt->description = "transferred from .nl file";
  //		//
  //		//				// count the number of entries
  //		//				if (d->kind & 4) // bit-wise mask to distinguish real from integer
  //		//						{
  //		//					conopt->type = "real";
  //		//					nOtherIdx = 0;
  //		//					for (int k = 0; k < n_con; k++) {
  //		//						if (d->u.r[k] != 0) {
  //		//							conopt->con[nOtherIdx] = new OtherConOption();
  //		//							conopt->con[nOtherIdx]->idx = k;
  //		//							conopt->con[nOtherIdx]->value = os_dtoa_format(
  //		//									d->u.r[k]);
  //		//							nOtherIdx++;
  //		//						}
  //		//					}
  //		//				} else             // here the suffix values are integer
  //		//				{
  //		//					conopt->type = "integer";
  //		//					nOtherIdx = 0;
  //		//					for (int k = 0; k < n_con; k++) {
  //		//						if (d->u.i[k] != 0) {
  //		//							conopt->con[nOtherIdx] = new OtherConOption();
  //		//							conopt->con[nOtherIdx]->idx = k;
  //		//							conopt->con[nOtherIdx]->value = os_dtoa_format(
  //		//									(double) d->u.i[k]);
  //		//							nOtherIdx++;
  //		//						}
  //		//					}
  //		//				}
  //		//
  //		//				conopt->numberOfCon = nOtherIdx;
  //		//
  //		//				if (found) {
  //		//					// here we just replace the <con> element and update the numberOfCon
  //		//					delete[] osoption->optimization->constraints->other[iopt]->con;
  //		//					osoption->optimization->constraints->other[iopt]->con =
  //		//							conopt->con;
  //		//					osoption->optimization->constraints->other[iopt]->numberOfCon =
  //		//							nOtherIdx;
  //		//					conopt->con = 0;
  //		//					osoption->optimization->constraints->other[iopt]->description =
  //		//							conopt->description;
  //		//				} else {
  //		//					if (!osoption->setAnOtherConstraintOption(conopt))
  //		//						throw ErrorClass(
  //		//								"OSnl2OS: Error transfering suffixes on constraints");
  //		//				}
  //		//
  //		//				delete conopt;
  //		//				conopt = 0;
  //		//			}
  //	}
  ////
  //	suffixType = ASL_Sufkind_obj;
  //	if ((asl->i.suffixes[suffixType] != 0)) {
  //		std::cout << "ya des suffixes" << __FILE__ << " l. " << __LINE__
  //				<< std::endl;
  ////				OtherObjectiveOption* objopt;
  //		//			for (d = asl->i.suffixes[suffixType]; d; d = d->next) {
  //		//#ifdef AMPLDEBUG
  //		//				std::cout << "Detected suffix " << d->sufname << "; kind = "
  //		//						<< d->kind << std::endl;
  //		//#endif
  //		//
  //		//				// allocate space
  //		//				objopt = new OtherObjectiveOption();
  //		//
  //		//				objopt->name = d->sufname;
  //		//				objopt->numberOfEnumerations = 0;
  //		//				objopt->obj = new OtherObjOption*[n_obj];
  //		//
  //		//				// check if the option was present in the OSoL file
  //		//				found = false;
  //		//				int iopt;
  //		//				for (iopt = 0; iopt < nOther; iopt++) {
  //		//					if (d->sufname == otherOptionNames[iopt]) {
  //		//						found = true;
  //		//						break;
  //		//					}
  //		//				}
  //		//
  //		//				// merge values by overwriting .nl file info
  //		//				if (found) {
  //		//					OtherObjectiveOption* otherOption;
  //		//					otherOption = osoption->getOtherObjectiveOption(iopt);
  //		//					for (int i = 0; i < otherOption->numberOfObj; i++) {
  //		//						if (d->kind & 4) // bit-wise mask to distinguish real from integer
  //		//								{
  //		//							d->u.r[otherOption->obj[i]->idx] = os_strtod(
  //		//									otherOption->obj[i]->value.c_str(), 0);
  //		//						} else
  //		//							d->u.i[otherOption->obj[i]->idx] = (int) os_strtod(
  //		//									otherOption->obj[i]->value.c_str(), 0)
  //		//									+ 0.1;
  //		//					}
  //		//					if (otherOption->description == "")
  //		//						objopt->description = "combined from osol and .nl data";
  //		//					else
  //		//						objopt->description = otherOption->description
  //		//								+ "; merged with .nl data";
  //		//				} else
  //		//					objopt->description = "transferred from .nl file";
  //		//
  //		//				// count the number of entries
  //		//				if (d->kind & 4) // bit-wise mask to distinguish real from integer
  //		//						{
  //		//					objopt->type = "real";
  //		//					nOtherIdx = 0;
  //		//					for (int k = 0; k < n_obj; k++) {
  //		//						if (d->u.r[k] != 0) {
  //		//							objopt->obj[nOtherIdx] = new OtherObjOption();
  //		//							objopt->obj[nOtherIdx]->idx = k;
  //		//							objopt->obj[nOtherIdx]->value = os_dtoa_format(
  //		//									d->u.r[k]);
  //		//							nOtherIdx++;
  //		//						}
  //		//					}
  //		//				} else             // here the suffix values are integer
  //		//				{
  //		//					objopt->type = "integer";
  //		//					nOtherIdx = 0;
  //		//					for (int k = 0; k < n_obj; k++) {
  //		//						if (d->u.i[k] != 0) {
  //		//							objopt->obj[nOtherIdx] = new OtherObjOption();
  //		//							objopt->obj[nOtherIdx]->idx = k;
  //		//							objopt->obj[nOtherIdx]->value = os_dtoa_format(
  //		//									(double) d->u.i[k]);
  //		//							nOtherIdx++;
  //		//						}
  //		//					}
  //		//				}
  //		//
  //		//				objopt->numberOfObj = nOtherIdx;
  //		//
  //		//				if (found) {
  //		//					// here we just replace the <var> element and update the numberOfVar
  //		//					delete[] osoption->optimization->objectives->other[iopt]->obj;
  //		//					osoption->optimization->objectives->other[iopt]->obj =
  //		//							objopt->obj;
  //		//					osoption->optimization->objectives->other[iopt]->numberOfObj =
  //		//							nOtherIdx;
  //		//					objopt->obj = 0;
  //		//					osoption->optimization->objectives->other[iopt]->description =
  //		//							objopt->description;
  //		//				} else {
  //		//					if (!osoption->setAnOtherObjectiveOption(objopt))
  //		//						throw ErrorClass(
  //		//								"OSnl2OS: Error transfering suffixes on objectives");
  //		//				}
  //		//
  //		//				delete objopt;
  //		//				objopt = 0;
  //		//			}
  //	}
  ////
  ////		//problem-indexed suffixes: is there ever more than one value in a .nl file?
  //	suffixType = ASL_Sufkind_prob;
  ////	if ((asl->i.suffixes[suffixType] != 0)) {
  ////
  ////		std::cout << "ya des suffixes" << __FILE__ << " l. " << __LINE__
  ////				<< std::endl;
  ////		std::string opttype, optvalue, optdesc;
  ////		optdesc = "transferred from .nl file";
  ////		for (d = asl->i.suffixes[suffixType]; d; d = d->next) {
  ////#ifdef AMPLDEBUG
  ////			std::cout << "Detected suffix " << d->sufname << "; kind = "
  ////					<< d->kind << std::endl;
  ////#endif
  ////
  ////			if (d->kind & 4) // bit-wise mask to distinguish real from integer
  ////					{
  ////				opttype = "real";
  ////				optvalue = os_dtoa_format(d->u.r[0]);
  ////			} else             // here the suffix values are integer
  ////			{
  ////				opttype = "integer";
  ////				optvalue = os_dtoa_format((double) d->u.i[0]);
  ////			}
  ////
  ////			if (!osoption->setAnotherSolverOption(d->sufname, optvalue, "", "",
  ////					opttype, optdesc))
  ////				throw ErrorClass(
  ////						"OSnl2OS: Error transfering problem-indexed suffixes");
  ////		}
  ////	}
  //
  //// process initial primal values
  //
  ////	std::cout << "Initial variable values:" << std::endl;
  ////	for (int i = 0; i < n_var; i++) {
  ////		std::cout << X0[i] << "    " << (int) havex0[i] << std::endl;
  ////	}
  //
  ////	// process and merge data contained in OSoL file
  ////	// .nl file has data in dense form, but we store in sparse form
  ////	// OSoL data supersede .nl file, so we turn redundant .nl file info OFF
  ////	if (osoption != 0 && osoption->optimization != 0
  ////			&& osoption->optimization->variables != 0
  ////			&& osoption->optimization->variables->initialVariableValues != 0
  ////			&& osoption->optimization->variables->initialVariableValues->numberOfVar
  ////					> 0) {
  ////		int n_prev =
  ////				osoption->optimization->variables->initialVariableValues->numberOfVar;
  ////		for (int i = 0; i < n_prev; i++)
  ////			havex0[osoption->optimization->variables->initialVariableValues->var[i]->idx] =
  ////					0;
  ////	}
  ////	std::cout << "Initial variable values, updated:" << std::endl;
  ////	for (int i = 0; i < n_var; i++) {
  ////		std::cout << X0[i] << "    " << (int) havex0[i] << std::endl;
  ////	}
  ////
  //// count the number of values (including those in the OSoL file)
  //	int n_x0 = 0;
  //	for (int i = 0; i < n_var; i++)
  //		if (havex0[i] != 0)
  //			n_x0++;
  //	p.initPoint().clear();
  //	if (n_x0 > 0) {
  //		// pull info out of ASL data structure
  //		n_x0 = 0;
  //		for (int i = 0; i < n_var; i++) {
  //			if (havex0[i] != 0) {
  //				p.initPoint()[i] = X0[i];
  //			}
  //		}
  //	}

  // process initial dual values

  //	std::cout << "Initial dual variable values:" << std::endl;
  //	for (int i = 0; i < n_con; i++) {
  //		std::cout << pi0[i] << "    " << (int) havepi0[i] << std::endl;
  //	}

  ////		// process and merge data contained in OSoL file
  ////		// .nl file has data in dense form, but we store in sparse form
  ////		// OSoL data supersede .nl file, so we turn redundant .nl file info OFF
  ////		if (osoption != 0 && osoption->optimization != 0
  ////				&& osoption->optimization->constraints != 0
  ////				&& osoption->optimization->constraints->initialDualValues
  ////						!= 0
  ////				&& osoption->optimization->constraints->initialDualValues->numberOfCon
  ////						> 0) {
  ////			int n_prev =
  ////					osoption->optimization->constraints->initialDualValues->numberOfCon;
  ////			for (int i = 0; i < n_prev; i++)
  ////				havepi0[osoption->optimization->constraints->initialDualValues->con[i]->idx] =
  ////						0;
  ////		}
  ////
  ////		std::cout << "Initial dual variable values, updated:" << std::endl;
  ////		for (int i = 0; i < n_con; i++) {
  ////			std::cout << pi0[i] << "    " << (int) havepi0[i] << std::endl;
  ////		}
  ////
  ////		// count the number of values (including those in the OSoL file)
  ////		int n_pi0 = 0;
  ////		for (int i = 0; i < n_con; i++)
  ////			if (havepi0[i] != 0)
  ////				n_pi0++;
  ////
  ////		if (n_pi0 > 0) {
  ////			if (osoption == 0) {
  ////				if (osol != "")
  ////					osoption = osolreader->readOSoL(osol);
  ////				else
  ////					osoption = new OSOption();
  ////			}
  ////
  ////			// allocate space
  ////			InitDualVarValue **pi_init = new InitDualVarValue*[n_pi0];
  ////
  ////			// pull info out of ASL data structure
  ////			n_pi0 = 0;
  ////			for (int i = 0; i < n_con; i++) {
  ////				if (havepi0[i] != 0) {
  ////					pi_init[n_pi0] = new InitDualVarValue();
  ////					pi_init[n_pi0]->idx = i;
  ////					pi_init[n_pi0]->lbDualValue = pi0[i];
  ////					pi_init[n_pi0]->ubDualValue = pi0[i];
  ////					n_pi0++;
  ////				}
  ////			}
  ////
  ////			// store into osoption
  ////			if (!osoption->setInitDualVarValuesSparse(n_pi0, pi_init,
  ////					ENUM_COMBINE_ARRAYS_merge))
  ////				throw ErrorClass(
  ////						"OSnl2OS: Error merging initial dual variable values");
  ////		}
  ////
  ////		//still to do
  ////		//initial basis status: .sstatus
  ////		/*
  ////		 0 - none - unknown
  ////		 1 - basic - basic
  ////		 2 - superbasic - superbasic
  ////		 3 - nonbasic at lower bound - atLower
  ////		 4 - nonbasic at upper bound - atUpper
  ////		 5 - nonbasic at equality bound - ?
  ////		 6 - nonbasic between bounds - isFree?
  ////		 */
  ////		//special ordered sets, branching weights, branching group weights
  ////		//initial objective values: .val
  ////	}    // end try
  ////
  ////	catch (const ErrorClass& eclass) {
  ////		// garbage collection etc.
  ////	}
  ////
  ////#ifdef AMPLDEBUG
  ////	OSiLWriter osilwriter;
  ////	std::cout << "WRITE THE INSTANCE" << std::endl;
  ////	osilwriter.m_bWhiteSpace = true;
  ////	std::cout << osilwriter.writeOSiL(osinstance) << std::endl;
  ////	std::cout << "DONE WRITE THE INSTANCE" << std::endl;
  ////
  ////	std::cout << osinstance->printModel() << std::endl;
  ////
  ////	OSoLWriter osolwriter;
  ////	std::cout << "WRITE THE OPTIONS" << std::endl;
  ////	//    osilwriter.m_bWhiteSpace = true;
  ////	std::cout << osolwriter.writeOSoL(osoption) << std::endl;
  ////	std::cout << "DONE WRITE THE OPTIONS" << std::endl;
  ////
  ////#endif
  //	for (auto & xp : nlExprs) {
  //		xp.check();
  //		//		std::cout << "xp : " << xp << std::endl;
  //		xp.free();
  //	}
  return true;
}
