 

#include "query/Operator.h"
#include "system/SystemCatalog.h"
#include "system/Exceptions.h"

using namespace std;

namespace scidb {

/**
 * @brief The operator: fuu().
 *
 * @par Synopsis:
 *   fuu( srcArray, probability [, seed] )
 *
 * @par Summary:
 *   Evaluates whether to include a cell in the result array by generating a random number and checks if it is less than probability.
 *
 * @par Input:
 *   - srcArray: a source array with srcAttrs and srcDims.
 *   - probability: the probability threshold, in [0..1]
 *   - an optional seed for the random number generator.
 *
 * @par Output array:
 *        <
 *   <br>   srcAttrs
 *   <br> >
 *   <br> [
 *   <br>   srcDims
 *   <br> ]
 *
 * @par Examples:
 *
 * @par Errors:
 *   n/a
 *
 * @par Notes:
 *   n/a
 *
 */
class LogicalDem: public LogicalOperator
{
public:
    LogicalDem(const string& logicalName, const std::string& alias):
        LogicalOperator(logicalName, alias)
    {
    	ADD_PARAM_INPUT()
        //ADD_PARAM_CONSTANT("double");
        //ADD_PARAM_VARIES()
    }
/*
	std::vector<std::shared_ptr<OperatorParamPlaceholder> > nextVaryParamPlaceholder(const std::vector< ArrayDesc> &schemas)
	{
	  std::vector<std::shared_ptr<OperatorParamPlaceholder> > res;
	  res.push_back(END_OF_VARIES_PARAMS());
	  res.push_back(PARAM_CONSTANT("int64"));
	  return res;
	}
*/

    ArrayDesc inferSchema(vector<ArrayDesc> schemas, std::shared_ptr< Query> query)
    {
        //assert(schemas.size() == 1);
        //return addEmptyTagAttribute(schemas[0]);
      
	// Build Output Array
	Attributes outputAttributes;
	outputAttributes.push_back( AttributeDesc(0, "t0", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(1, "t1", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(2, "t2", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(3, "t3", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(4, "t4", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(5, "t5", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(6, "t6", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(7, "t7", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(8, "t8", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(9, "t9", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(10, "t10", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(11, "t11", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(12, "t12", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(13, "t13", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(14, "t14", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(15, "t15", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(16, "t16", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(17, "t17", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(18, "t18", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(19, "t19", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(20, "t20", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(21, "t21", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(22, "t22", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(23, "t23", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(24, "t24", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(25, "t25", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(26, "t26", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(27, "t27", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(28, "t28", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(29, "t29", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(30, "t30", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(31, "t31", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(32, "t32", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(33, "t33", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(34, "t34", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(35, "t35", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(36, "t36", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(37, "t37", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(38, "t38", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
	outputAttributes.push_back( AttributeDesc(39, "t39", TID_DOUBLE, AttributeDesc::IS_NULLABLE, 0));
      
        return ArrayDesc("dem", outputAttributes,
                         schemas[0].getDimensions(),
                         schemas[0].getDistribution(),
                         schemas[0].getResidency());
      
    }
};

REGISTER_LOGICAL_OPERATOR_FACTORY(LogicalDem, "dem");


} //namespace