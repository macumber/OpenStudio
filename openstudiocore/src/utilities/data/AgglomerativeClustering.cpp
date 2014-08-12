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

  ClusteringResult::ClusteringResult(const std::string& name, double score, const std::vector<Cluster>& clusters)
    : m_name(name), m_score(score), m_clusters(clusters)
  {}

  std::string ClusteringResult::name() const
  {
    return m_name;
  }

  double ClusteringResult::score() const
  {
    return m_score;
  }

  std::vector<Cluster> ClusteringResult::clusters() const
  {
    return m_clusters;
  }

  AgglomerativeClusterer::AgglomerativeClusterer(const TimeSeries& timeSeries)
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
    } catch (const std::exception& ){
      LOG_AND_THROW("Cannot partition TimeSeries into Vectors of equal length");
    }
  }

  struct BestMerge
  {
    unsigned i;
    unsigned j;
    Cluster cluster;
  };

  bool AgglomerativeClusterer::solve()
  {
    // n is number of vectors, p is number of clusters in result
    unsigned n = m_vectors.size();
    unsigned p = n;
    double ss0 = 0;

    // start with each vector in its own cluster
    std::vector<Cluster> clusters;
    clusters.reserve(m_vectors.size());
    for (const Vector& vector : m_vectors){
      Cluster cluster(vector);
      ss0 += cluster.sumOfSquares();
      clusters.push_back(cluster);
    }
    ClusteringResult initialClustering("Initial Clustering", ss0, clusters);
    m_clusteringResults.push_back(initialClustering);

    // now reduce p by 1 until you get to a single cluster
    // matrix ss is p x p, where ss[i,j] is the ss from merging clusters[i] with clusters[j]
    Matrix ss(p, p, std::numeric_limits<double>::max());
    for (unsigned i = 0; i < p; ++i){
    }

      struct BestMerge
      {
        unsigned i;
        unsigned j;
        double ss;
      };
    


    return true;
  }

  void AgglomerativeClusterer::clear()
  {
    m_clusterResults.clear();
  }

  std::vector<ClusteringResult> AgglomerativeClusterer::clusterResults() const
  {
    std::vector<ClusteringResult> result(m_clusterResults);
    result.insert(result.end(), m_specialClusterResults.begin(), m_specialClusterResults.end());
    return result;
  }

} // openstudio
