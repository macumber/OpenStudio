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

  ClusterData::ClusterData(const std::vector<Vector>& vectors)
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
      }
      else {
        n = vector.size();
        if (n.get() == 0){
          LOG_AND_THROW("Vectors must have length greater than 0");
        }
      }
    }
  }

  unsigned ClusterData::numVectors() const
  {
    return m_vectors.size();
  }

  std::string ClusterData::units() const
  {
    return m_units;
  }

  const std::vector<Vector>& ClusterData::vectors() const
  {
    return m_vectors;
  }

  const std::vector<Time>& ClusterData::times() const
  {
    return m_times;
  }

  const std::vector<Date>& ClusterData::dates() const
  {
    return m_dates;
  }

  void ClusterData::setUnits(const std::string& units)
  {
    m_units = units;
  }

  bool ClusterData::setTimes(const std::vector<Time>& times)
  {
    OS_ASSERT(!m_vectors.empty());

    bool result = false;
    if (times.size() == m_vectors[0].size()){
      m_times = times;
      result = true;
    }

    return result;
  }

  bool ClusterData::setDates(const std::vector<Date>& dates)
  {
    bool result = false;
    if (dates.size() == m_vectors.size()){
      m_dates = dates;
      result = true;
    }

    return result;
  }

  //std::string ClusterData::indicesToKey(const std::vector<unsigned>& indices) const
  //{
  //  std::stringstream ss;
  //  for (unsigned i : indices){
  //    ss << i << ",";
  //  }
  //  return ss.str();
  //}

  boost::optional<Cluster> ClusterData::getCachedCluster(const std::vector<unsigned>& indices)
  {
    //std::map<std::string, Cluster>::const_iterator it = m_cachedClusters.find(indicesToKey(indices));
    std::map<std::vector<unsigned>, Cluster>::const_iterator it = m_cachedClusters.find(indices);
    if (it != m_cachedClusters.end()){
      return it->second;
    }
    return boost::none;
  }

  void ClusterData::addCachedCluster(const std::vector<unsigned>& indices, const Cluster& cluster)
  {
   // m_cachedClusters.insert(std::make_pair<std::string, Cluster>(indicesToKey(indices), Cluster(cluster)));
    m_cachedClusters.insert(std::make_pair<std::vector<unsigned>, Cluster>(std::vector<unsigned>(indices), Cluster(cluster)));
  }

  Cluster::Cluster(std::shared_ptr<ClusterData> clusterData, const std::vector<unsigned>& indices)
    : m_clusterData(clusterData), m_indices(indices)
  {
    OS_ASSERT(clusterData);

    const std::vector<Vector>& allVectors = clusterData->vectors();
    unsigned n = allVectors.size();
    unsigned m = indices.size();
    for (unsigned i : indices){
      if (i > (n - 1)){
        LOG_AND_THROW("Index out of range");
      }
      Vector vector = allVectors[i];      
      m_vectors.push_back(vector);
      Vector scaledVector = (1.0 / m)*vector;
      if (m_meanVector.empty()){
        m_meanVector = scaledVector;
      }else{
        m_meanVector += scaledVector;
      }
    }

    m_sumOfSquares = 0;
    for (const Vector& vector : m_vectors){
      Vector error = vector - m_meanVector;
      m_sumOfSquares += pow(boost::numeric::ublas::norm_2(error), 2.0);
    }

    m_clusterData->addCachedCluster(indices, *this);
  }

  const std::vector<unsigned>& Cluster::indices() const
  {
    return m_indices;
  }

  const std::vector<Vector>& Cluster::vectors() const
  {
    return m_vectors;
  }

  unsigned Cluster::numVectors() const
  {
    return m_vectors.size();
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
    std::vector<unsigned> newIndices(m_indices);
    std::vector<unsigned> otherIndices(other.indices());
    newIndices.insert(newIndices.end(), otherIndices.begin(), otherIndices.end());
    std::sort(newIndices.begin(), newIndices.end());

    boost::optional<Cluster> result = m_clusterData->getCachedCluster(newIndices);
    if (result){
      return *result;
    }

    return Cluster(m_clusterData, newIndices);
  }

  ClusteringResult::ClusteringResult(const std::string& name, const std::vector<Cluster>& clusters, double singleClusterSumOfSquares)
    : m_name(name), m_clusters(clusters), m_singleClusterSumOfSquares(singleClusterSumOfSquares), m_sumOfSquares(0)
  {
    unsigned n = 0; // total number of vectors
    unsigned p = clusters.size(); // number of clusters
    for (const Cluster& cluster : m_clusters){
      n += cluster.numVectors();
      m_sumOfSquares += cluster.sumOfSquares();
    }

    m_rSquared = (m_singleClusterSumOfSquares - m_sumOfSquares) / m_singleClusterSumOfSquares;
    if (n == p){
      m_rSquaredAdjusted = 0;
    }else{
      m_rSquaredAdjusted = m_rSquared - (1 - m_rSquared)*(p / (n - p));
    }
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
    OS_ASSERT(bestI);
    OS_ASSERT(bestJ);

    std::stringstream ss;
    ss << "Clustering " << n - 1;
    std::string name = ss.str();

    std::vector<Cluster> newClusters;
    for (unsigned i = 0; i < n; ++i){
      if (i == bestI.get()){
        newClusters.push_back(*bestCluster);
      }else if (i == bestJ.get()){
        // merged with bestI
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

    std::vector<std::vector<double> > stdVectors;
    stdVectors.reserve(n);

    std::vector<Time> times;
    std::vector<Date> dates;
    unsigned timeIndex = 0;

    for (unsigned i = 0; i < n; ++i){
      if (!dates.empty()){
        if (dateTimes[i].date() != dates.back() && dateTimes[i].time().totalHours() > 0.0){
          // new date
          timeIndex = 0;
          if (timeIndex > (times.size()-1) || times[timeIndex] != dateTimes[i].time()){
            LOG_AND_THROW("Time indices do not match.");
          }
          dates.push_back(dateTimes[i].date());
          stdVectors.push_back(std::vector<double>());
          stdVectors.back().reserve(n);
        } else if (dates.size() == 1){
          // still on first day
          times.push_back(dateTimes[i].time());
        }else{
          if (timeIndex > (times.size() - 1) || times[timeIndex] != dateTimes[i].time()){
            LOG_AND_THROW("Time indices do not match.");
          }
        }
      }else{
        // first point processed
        times.push_back(dateTimes[i].time());
        dates.push_back(dateTimes[i].date());
        stdVectors.push_back(std::vector<double>());
        stdVectors.back().reserve(n);
      }
      
      ++timeIndex;
      stdVectors.back().push_back(values[i]);
    }

    n = stdVectors.size();

    std::vector<Vector> vectors;
    vectors.reserve(n);
    for (const std::vector<double>& stdVector : stdVectors){
      vectors.push_back(createVector(stdVector));
    }

    std::vector<unsigned> allIndices;
    allIndices.reserve(n);
    for (unsigned i = 0; i < n; ++i){
      allIndices.push_back(i);
    }

    // this can throw
    m_clusterData = std::shared_ptr<ClusterData>(new ClusterData(vectors));
    m_clusterData->setUnits(timeSeries.units());
    m_clusterData->setDates(dates);
    m_clusterData->setTimes(times);

    // this can throw
    Cluster test(m_clusterData, allIndices);
    m_singleClusterSumOfSquares = test.sumOfSquares();

  }

  bool AgglomerativeClusterer::solve()
  {
    unsigned n = m_clusterData->numVectors();

    // start with each vector in its own cluster
    std::vector<Cluster> clusters;
    clusters.reserve(n);
    for (unsigned i = 0; i < n; ++i){
      std::vector<unsigned> indices(1, i);
      Cluster cluster(m_clusterData, indices);
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
