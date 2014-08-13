/**********************************************************************
*  Copyright (c) 2008-2014, Alliance for Sustainable Energy.  
*  All rights reserved.
*  
*  This library is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Lesser General Public
*  License as published by the Free Software Foundation; either
*  version 2.1 of the License, or (at your option) any later version.
*  
*  This library is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*  Lesser General Public License for more details.
*  
*  You should have received a copy of the GNU Lesser General Public
*  License along with this library; if not, write to the Free Software
*  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
**********************************************************************/

#include <gtest/gtest.h>
#include "DataFixture.hpp"

#include "../AgglomerativeClustering.hpp"
#include "../TimeSeries.hpp"
#include "../Vector.hpp"
#include "../../time/Date.hpp"
#include "../../time/Time.hpp"

using namespace std;
using namespace boost;
using namespace openstudio;

TEST_F(DataFixture, AgglomerativeClustering_TimeSeriesConstructor)
{

  Vector values = randVector(-1.0, 1.0, 8760);
  TimeSeries timeSeries(Date(1, 1), Time(0, 1), values, "W");

  AgglomerativeClusterer clusterer(timeSeries);
  clusterer.solve();

  std::vector<ClusteringResult> clusteringResults = clusterer.clusteringResults();
  for (const ClusteringResult& clusteringResult : clusteringResults){
    std::cout << clusteringResult.name() << " " << clusteringResult.sumOfSquares() << " " << clusteringResult.rSquared() << " " << clusteringResult.rSquaredAdjusted() << std::endl;
  }

}