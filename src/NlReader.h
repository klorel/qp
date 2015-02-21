/* $Id: OSnl2OS.h 4249 2011-08-11 01:08:14Z Gassmann $ */
/** @file OSnl2osol.h
 * 
 * @author  Horand Gassmann, Jun Ma, Kipp Martin
 *
 * \remarks
 * Copyright (C) 2011, Horand Gassmann, Jun Ma, Kipp Martin,
 * Dalhousie University, and the University of Chicago.
 * All Rights Reserved.
 * This software is licensed under the Eclipse Public License. 
 * Please see the accompanying LICENSE file in root directory for terms.
 * 
 */

#ifndef NLREADER_2_H
#define NLREADER_2_H

#include <string>
#include <vector>
#include <map>
#include <cstdint>

#include "Problem.h"

struct ograd;
struct cgrad;
struct ASL;
struct expr;

enum OS_AMPL_SUFFIX_TYPE {
	OS_AMPL_SUFFIX_TYPE_integer, OS_AMPL_SUFFIX_TYPE_double
};

enum OS_AMPL_SUFFIX_SCOPE {
	OS_AMPL_SUFFIX_SCOPE_variables,
	OS_AMPL_SUFFIX_SCOPE_constraints,
	OS_AMPL_SUFFIX_SCOPE_objectives,
	OS_AMPL_SUFFIX_SCOPE_problems
};

enum OS_AMPL_SUFFIX_DIRECTION {
	OS_AMPL_SUFFIX_DIRECTION_toSolver,
	OS_AMPL_SUFFIX_DIRECTION_toAMPL,
	OS_AMPL_SUFFIX_DIRECTION_both,
	OS_AMPL_SUFFIX_DIRECTION_local
};

struct OS_AMPL_SUFFIX {
	std::string name;
	OS_AMPL_SUFFIX_TYPE type;
	OS_AMPL_SUFFIX_SCOPE scope;
	OS_AMPL_SUFFIX_DIRECTION direction;
};

/*
 void OS_addAmplSuffix(std::string name,
 OS_AMPL_SUFFIX_TYPE type,
 OS_AMPL_SUFFIX_SCOPE scope,
 OS_AMPL_SUFFIX_DIRECTION direction)
 {
 OS_AMPL_SUFFIX temp;
 temp.name = name;
 temp.type = type;
 temp.scope = scope;
 temp.direction = direction;
 OS_amplSuffixTable.push_back(temp);
 };
 */

class NlReader {
public:
	/** the OSnl2OS class constructor */
	NlReader(std::string nlfilename, std::string osol);

	/** the OSnl2OS class destructor */
	~NlReader();

	/**
	 * create an OSInstance and OS option representation from the AMPL nl content
	 * Since some of the information in the nl file
	 * (such as initial values, basis information, branching priorities)
	 * cannot be stored into an OSInstance and must be stored in an
	 * OSOption object instead, all other options must be provided
	 * to this method as well.
	 *
	 * @param osol is a string containing the option information
	 * (e.g., as read from an osol file). 
	 * @return whether the OSOption object is created successfully.
	 * @remark The osol string must be prepared prior to calling createOSObjects!
	 */
	bool createObjects(Problem &); //(std::string osol);
	bool createVariables(Problem &);
	bool createConstraints(Problem &);
	bool createObjectif(Problem &);

	/**
	 * parse an nl tree structure holding a nonlinear expression
	 *
	 * @return the AMPL nonlinear structure as an OSnLNode.
	 */
//    OSnLNode* walkTree(expr *e);
	Function walkTree(Problem & p, expr *e);

	static std::map<size_t, std::string> OppName;
	static void BuildOppName();

	/** osol is a string containing the content of the OS option file
	 *  (it may be empty if no option file was provided). 
	 *  If osoption is 0, the option information is found in osol.
	 */
	std::string osol;

	std::vector<std::string> op_type;
	std::vector<double> operand;
	int numkount;
	bool _onoff;

private:

	/** og is a pointer to the AMPL data structure holding the
	 * objective function coefficients
	 */
	ograd *og;

	/** Pointers to AMPL data structures.
	 * cw is loaded in column-wise format.
	 * rw is loaded in row-wise format.
	 * asl is for conveniently switching.
	 */
	ASL *asl;

	void getRW();
	void getCW();
	void free();

	/** stub is the name of the file with the nl instance
	 */
	std::string stub;
	bool _isCW;
};
//end of  OSnl2OS

#endif
