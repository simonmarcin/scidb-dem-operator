/*
**
* BEGIN_COPYRIGHT
*
* Copyright (C) 2008-2015 SciDB, Inc.
* All Rights Reserved.
*
* SciDB is free software: you can redistribute it and/or modify
* it under the terms of the AFFERO GNU General Public License as published by
* the Free Software Foundation.
*
* SciDB is distributed "AS-IS" AND WITHOUT ANY WARRANTY OF ANY KIND,
* INCLUDING ANY IMPLIED WARRANTY OF MERCHANTABILITY,
* NON-INFRINGEMENT, OR FITNESS FOR A PARTICULAR PURPOSE. See
* the AFFERO GNU General Public License for the complete license terms.
*
* You should have received a copy of the AFFERO GNU General Public License
* along with SciDB.  If not, see <http://www.gnu.org/licenses/agpl-3.0.html>
*
* END_COPYRIGHT
*/

/*
 * @file PhysicalDem.cpp
 * @author simon.marcin@fhnw.ch
 *
 * @brief Structure to loop through all Chunks and getAttributes
 * 	  and assign them to a AlgorithmObject.
 */

#include "query/Operator.h"
#include <iostream>
#include <fstream>
#include <Firdem2.h>

using namespace std;
using namespace boost;

namespace scidb
{

  
class PhysicalDem: public PhysicalOperator
{
public:
    PhysicalDem(const std::string& logicalName, const std::string& physicalName, const Parameters& parameters, const ArrayDesc& schema):
	    PhysicalOperator(logicalName, physicalName, parameters, schema)
	{
	}

    std::shared_ptr<Array> execute(std::vector<std::shared_ptr<Array> >& inputArrays, std::shared_ptr<Query> query)
	{

	  
	  //Init firdem
	  Firdem2 firdem;
	  boost::numeric::ublas::vector<double> coffs(40);
	  boost::numeric::ublas::vector<int> data(6);
	  int itcount;

	  //Init of Array and Chunk Iterators
	  size_t nAttrs = inputArrays[0]->getArrayDesc().getAttributes().size();
	  vector<std::shared_ptr<ConstArrayIterator> > srcArrayIterators (nAttrs);
	  vector<std::shared_ptr<ConstChunkIterator> > srcChunkIterators (nAttrs);
	  
	  //Debug files
	  //InstanceID instanceId = query->getInstanceID();
	  //ofstream myfile;
	  //std::string file ="/tmp/dem_log.";
	  //file += boo st::lexical_cast<std::string>(instanceId);
	  //myfile.open(file);
	  //myfile << "Debug Log.\n";

	  
	  //measure time
	  //time_t now = time(0);
	  //char* dt = ctime(&now);
	  //myfile << dt << " - Physical Operator started\n";
	
	//Init of Array and Chunk Iterators for Output Array
	std::shared_ptr<Array> outputArray(new MemArray(_schema, query));
	size_t outNAttrs = outputArray->getArrayDesc().getAttributes().size();
	vector<std::shared_ptr<ArrayIterator> > outArrayIterators (outNAttrs);
	vector<std::shared_ptr<ChunkIterator> > outChunkIterators (outNAttrs);
	
	//Get ArrayIterator for each attribute
	for (size_t i=0;i<nAttrs;i++){
	  srcArrayIterators[i] = inputArrays[0]->getConstIterator(i);
	}
	
	//Get ArrayIterator for each attribute - OutputArray
	for (size_t i=0;i<outNAttrs;i++){
	  outArrayIterators[i] = outputArray->getIterator(i);
	}
	
	//Process all Chunks
	while (!srcArrayIterators[0]->end())
	{
	  
	  //Get ChunkIterator for each Attribute
	  for(size_t i=0;i<nAttrs;i++)
	  {
	    ConstChunk const& inputChunk = srcArrayIterators[i]->getChunk();
	    srcChunkIterators[i] = inputChunk.getConstIterator();
	    //srcChunkIterators[i] = srcArrayIterators[i]->getChunk().getConstIterator(0);
	  }
	  
	  //Same here for the OutputArray
	  for(size_t i=0;i<outNAttrs;i++)
	  {
	    Coordinates position = srcChunkIterators[0]->getPosition();
	    outChunkIterators[i] = outArrayIterators[i]->newChunk(position).getIterator(query,
                                                                     ChunkIterator::SEQUENTIAL_WRITE);
	    outChunkIterators[i]->setPosition(position);
	  }
	  
	    //Dem through the Chunks
	    while (!srcChunkIterators[0]->end())
	    {
	      
	      //For each attribute
	      for (size_t i=0;i<nAttrs-1;i++)
	      {
		
		Value const& inputValue = srcChunkIterators[i]->getItem();
		data(i) = inputValue.getInt32();
		Value val = inputValue;
		
		double tmp = val.getDouble();
		tmp=tmp+2.2;
		val.setDouble(tmp);
		
		
		++(*srcChunkIterators[i]);
			
               }
	       
	       firdem.calculateDem(data,coffs,itcount);
               
	       
	      //Write all Attributes to the Output Chunk
	      for (size_t i=0;i<outNAttrs;i++)
	      {
		Value value;
		value.setDouble(coffs(i));
		outChunkIterators[i]->writeItem(value);
		++(*outChunkIterators[i]);
	      }
                
            }
            
            //finish output chunk
            for (size_t i=0;i<outNAttrs;i++)
	    {
	      outChunkIterators[i]->flush();
	      ++(*outArrayIterators[i]);
	    }
            
            
	  //Go to the next Chunk
	  for(size_t i=0;i<nAttrs;i++)
	  {
	    ++(*srcArrayIterators[i]);
	  }
	    
	}
	

	//now = time(0);
	//dt = ctime(&now);
	//myfile << dt << " - Physical Operator finished\n";
	//myfile.close();
	
        //return Output Array
        return outputArray;
        
	}
};

REGISTER_PHYSICAL_OPERATOR_FACTORY(PhysicalDem, "dem", "dem_impl");

} //namespace
