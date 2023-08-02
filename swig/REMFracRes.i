 //  ligne de commande : swig -Wall -c++ -java -outdir java/src/fracresSim/ -package fracresSim -o swig/fracresSim_wrap.cxx swig/fracresSim.i

%module REMFracRes

/* -----------------------------------------------------------------------------
 *
 *
 * REMFracRes.i
 *
 * ----------------------------------------------------------------------------- */
/* -----------------------------------------------------------------------------
 * See the LICENSE file for information on copyright, usage and redistribution
 * of SWIG, and the README file for authors - http://www.swig.org/release.html.
 *
 * REMFracRes.i
 *
 * ----------------------------------------------------------------------------- */

%{
	#include "/home/pascal/git/REMFracRes/include/remfracres/DistributionPower.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/fracResFileKeyWord.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResUnit.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResDataTypeEnum.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResCategoryEnum.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResClustersFaultTypeEnum.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResTypeEnum.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResDistributionArrayTypeEnum.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResDistributionEnum.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResGeometryEnum.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResGeometryCorrelationTypeEnum.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResBeddingAperture.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResBeddingProperties.h"	
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResFracture.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/BoxAABBEvalDistance.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/CellId.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FacetId.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/fractureset.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FaultArraySet.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResHypothesis.h"
    #include "/home/pascal/git/REMFracRes/include/remfracres/FracResScenario.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/FracResStage.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/RemFracResFileManagement.h"
	#include "/home/pascal/git/REMFracRes/include/remfracres/RemFracResSim.h"
	
	using namespace geode;
	using namespace boost;
	using namespace std;
	
%}

%include "std_vector.i"
%include "std_map.i"
%include "std_string.i"

// Instantiate templates used by example
namespace std {

%template(BoolVector) vector< bool >;
%template(Bool2DVector) vector< vector < bool > >;
%template(IntVector) vector< int >;
%template(Int2DVector) vector< vector < int > >;
%template(DoubleVector) vector < double >;
%template(FloatVector) vector < float >;
%template(Double2DVector) vector< vector < double > >;
%template(Double3DVector) vector< vector < vector < double > > >;
%template(StringVector) vector< string >;
%template(String2DVector) vector< vector< string > >;
%template(String3DVector) vector<vector <vector< string > > >;

}
%include "/home/pascal/git/REMFracRes/include/remfracres/DistributionPower.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/fracResFileKeyWord.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResUnit.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResDataTypeEnum.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResCategoryEnum.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResClustersFaultTypeEnum.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResTypeEnum.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResDistributionArrayTypeEnum.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResDistributionEnum.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResGeometryEnum.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResGeometryCorrelationTypeEnum.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResBeddingAperture.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResBeddingProperties.h"	
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResFracture.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/BoxAABBEvalDistance.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/CellId.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FacetId.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/fractureset.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FaultArraySet.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResHypothesis.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResScenario.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/FracResStage.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/RemFracResFileManagement.h"
%include "/home/pascal/git/REMFracRes/include/remfracres/RemFracResSim.h"



