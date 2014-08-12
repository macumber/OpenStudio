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

#include "AgglomerativeClustering.hpp"
#include "TimeSeries.hpp"
#include "../core/Assert.hpp"

#include <exception>

using namespace std;
using namespace boost;

namespace openstudio{

  Cluster::Cluster(const std::vector<Vector>& vectors)
    : m_vectors(vectors)
  {
    if (vectors.empty()){
      LOG_AND_THROW("Cannot create cluster with no Vectors");
    }

    boost::optional<unsigned> n;
    for (const Vector& vector : vectors){
      if (n){
        if (n.get() != vector.size()){
          LOG_AND_THROW("Vectors must be equal length");
        }
        m_meanVector += (vector / n.get());
      } else {
        n = vector.size();
        if (n.get() == 0){
          LOG_AND_THROW("Vectors must have length greater than 0");
        }else{
          m_meanVector = (vector / n.get());
        }
      }
    }

    m_sumOfSquares = 0;
    for (const Vector& vector : vectors){
      Vector error = vector - m_meanVector;
      m_sumOfSquares += pow(boost::numeric::ublas::norm_2(error), 2.0);
    }
    
  }

  std::vector<Vector> Cluster::vectors() const
  {
    return m_vectors;
  }

  Vector Cluster::meanVector() const
  {
    return m_meanVector;
  }

  double Cluster::sumOfSquares() const
  {
    return m_sumOfSquares;
  }

  Cluster Cluster::merge(const Cluster& other) const
  {
    std::vector<Vector> newVector(other.vectors());
    newVector.insert(newVector.end(), m_vectors.begin(), m_vectors.end());
    return Cluster(newVector);
  }

  ClusteringResult::ClusteringResult(const std::string& name, const std::vector<Cluster>& clusters, double singleClusterSumOfSquares)
    : m_name(name), m_clusters(clusters), m_singleClusterSumOfSquares(singleClusterSumOfSquares), m_sumOfSquares(0)
  {
    unsigned n = 0; // total number of vectors
    unsigned p = clusters.size(); // number of clusters
    for (const Cluster& cluster : m_clusters){
      n += cluster.vectors().size();
      m_sumOfSquares += cluster.sumOfSquares();
    }

    m_rSquared = (m_singleClusterSumOfSquares - m_sumOfSquares) / m_singleClusterSumOfSquares;
    m_rSquaredAdjusted = m_rSquared - (1 - m_rSquared)*(p / (n - p));
  }

  std::string ClusteringResult::name() const
  {
    return m_name;
  }

  double ClusteringResult::sumOfSquares() const
  {
    return m_sumOfSquares;
  }

  double ClusteringResult::rSquared() const
  {
    return m_rSquared;
  }

  double ClusteringResult::rSquaredAdjusted() const
  {
    return m_rSquaredAdjusted;
  }

  std::vector<Cluster> ClusteringResult::clusters() const
  {
    return m_clusters;
  }

  boost::optional<ClusteringResult> ClusteringResult::nextClusteringResult() const
  {
    boost::optional<unsigned> bestI;
    boost::optional<unsigned> bestJ;
    boost::optional<Cluster> bestCluster;

    unsigned n = m_clusters.size();
    for (unsigned i = 0; i < n; ++i){
      for (unsigned j = i + 1; j < n; ++j){
        Cluster temp = m_clusters[i].merge(m_clusters[j]);
        if (!bestCluster){
          bestI = i;
          bestJ = j;
          bestCluster = temp;
        }else if (temp.sumOfSquares() < bestCluster->sumOfSquares()){
          bestI = i;
          bestJ = j;
          bestCluster = temp;
        }
      }
    }

    if (!bestCluster){
      return boost::none;
    }

    std::stringstream ss;
    ss << "Clustering " << n - 1;
    std::string name = ss.str();
      
    std::vector<Cluster> newClusters;
    for (unsigned i = 0; i < n; ++i){
      if (i == bestI.get()){
        newClusters.push_back(*bestCluster);
      }else{
        newClusters.push_back(m_clusters[i]);
      }
    }
    OS_ASSERT(newClusters.size() == n - 1);

    return ClusteringResult(name, newClusters, m_singleClusterSumOfSquares);
  }

  AgglomerativeClusterer::AgglomerativeClusterer(const TimeSeries& timeSeries)
    : m_singleClusterSumOfSquares(0)
  {
    DateTimeVector dateTimes = timeSeries.dateTimes();
    openstudio::Vector values = timeSeries.values();
    unsigned n = values.size();
    OS_ASSERT(dateTimes.size() == n);

    boost::optional<Date> currentDate;
    std::vector<std::vector<double> > stdVectors;
    stdVectors.reserve(n);
    for (unsigned i = 0; i < n; ++i){
      if (currentDate){
        if (dateTimes[i].date() != *currentDate){
          stdVectors.push_back(std::vector<double>());
          stdVectors.back().reserve(n);
        }
      }else{
        currentDate = dateTimes[i].date();
        stdVectors.push_back(std::vector<double>());
        stdVectors.back().reserve(n);
      }
      
      stdVectors.back().push_back(values[i]);
    }

    m_vectors;
    m_vectors.reserve(n);
    for (const std::vector<double>& stdVector : stdVectors){
      m_vectors.push_back(createVector(stdVector));
    }

    try{
      // compute a test cluster to ensure these Vectors are all the same size
      Cluster test(m_vectors);
      m_singleClusterSumOfSquares = test.sumOfSquares();
    } catch (const std::exception& ){
      LOG_AND_THROW("Cannot partition TimeSeries into Vectors of equal length");
    }
  }

  bool AgglomerativeClusterer::solve()
  {
    // start with each vector in its own cluster
    std::vector<Cluster> clusters;
    clusters.reserve(m_vectors.size());
    for (const Vector& vector : m_vectors){
      std::vector<Vector> tmp(1, vector);
      Cluster cluster(tmp);
      clusters.push_back(cluster);
    }
    
    boost::optional<ClusteringResult> clusteringResult = ClusteringResult("Initial Clustering", clusters, m_singleClusterSumOfSquares);
    while (clusteringResult) {
      m_clusteringResults.push_back(*clusteringResult);
      clusteringResult = clusteringResult->nextClusteringResult();
    }

    return true;
  }

  void AgglomerativeClusterer::clear()
  {
    m_clusteringResults.clear();
  }

  std::vector<ClusteringResult> AgglomerativeClusterer::clusteringResults() const
  {
    return m_clusteringResults;
  }

  /// returns the special ClusteringResults
  std::vector<ClusteringResult> AgglomerativeClusterer::specialClusteringResults() const
  {
    return m_specialClusteringResults;
  }

} // openstudio
