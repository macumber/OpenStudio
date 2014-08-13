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

#ifndef UTILITIES_DATA_AGGLOMERATIVECLUSTERING_HPP
#define UTILITIES_DATA_AGGLOMERATIVECLUSTERING_HPP

#include "../UtilitiesAPI.hpp"

#include "Vector.hpp"
#include "../core/Logger.hpp"

#include <boost/optional.hpp>

#include <vector>

namespace openstudio{

  class Time;
  class Date;
  class TimeSeries;
  class AgglomerativeClusterer;
  class ClusteringResult;
  class Cluster;

  /** ClusterData holds the shared data between several clusters.
  **/
  class UTILITIES_API ClusterData
  {
  public:
    /** @name Constructors */
    //@{

    /// constructor from Vectors, throws if all Vectors are not the same length or if length is 0
    ClusterData(const std::vector<Vector>& vectors);

    //@}
    /** @name Getters */
    //@{

    /// returns the number of Vectors in this ClusterData
    unsigned numVectors() const;

    /// returns the units associated with each value
    std::string units() const;

    /// returns the Vectors in this ClusterData
    const std::vector<Vector>& vectors() const;

    /// returns the Times in this ClusterData
    const std::vector<Time>& times() const;

    /// returns the Dates in this ClusterData
    const std::vector<Date>& dates() const;

    //@}
    /** @name Setters */
    //@{

    /// sets the units associated with each value
    void setUnits(const std::string& units);

    /// sets the Times associated with each Vector, returns false if not the same length as all Vectors
    bool setTimes(const std::vector<Time>& times);

    /// sets the Dates associated with each Vector, returns false if not the same length as number of Vectors
    bool setDates(const std::vector<Date>& dates);

    //@}

    boost::optional<Cluster> getCachedCluster(const std::vector<unsigned>& indices);

    void addCachedCluster(const std::vector<unsigned>& indices, const Cluster& cluster);

  private:

    REGISTER_LOGGER("utilities.Cluster");

    std::string m_units;
    std::vector<Vector> m_vectors;
    std::vector<Time> m_times;
    std::vector<Date> m_dates;

    //std::string indicesToKey(const std::vector<unsigned>& indices) const;

    // cache previous clusters
    //std::map<std::string, Cluster> m_cachedClusters;
    std::map<std::vector<unsigned>, Cluster> m_cachedClusters;
  };


  /** A Cluster is a vector of Vectors that have been combined together.
  **/
  class UTILITIES_API Cluster
  {

    /// returns the indices in this cluster
    const std::vector<unsigned>& indices() const;

    /// returns the Vectors in this cluster
    const std::vector<Vector>& vectors() const;

    /// returns the number of Vectors in this cluster
    unsigned numVectors() const;

    /// computes the mean vector, $\hat x_{j} = \frac{1}{n} \sum_{i=0}^n x_{j}
    Vector meanVector() const;

    /// computes the within cluster sum of square error, $ss = \frac{1}{n} \sum_{i=0}^n \sum_{j=0}^m (x_{j}-\hat x_{j})^2
    double sumOfSquares() const;

    /// merges this cluster with another one, throws if all vectors are not the same length
    Cluster merge(const Cluster& other) const;

  private:

    REGISTER_LOGGER("utilities.Cluster");

    friend class AgglomerativeClusterer;
    friend class ClusteringResult;
    
    /// constructor from ClusterData, throws if any index is out of range
    Cluster(std::shared_ptr<ClusterData> clusterData, const std::vector<unsigned>& indices);

    std::shared_ptr<ClusterData> m_clusterData;
    std::vector<unsigned> m_indices;
   
    // cached
    std::vector<Vector> m_vectors;
    Vector m_meanVector;
    double m_sumOfSquares;
  };

  /** ClusteringResult holds a vector of Clusters which partition the entire space.
   **/
  class UTILITIES_API ClusteringResult
  {
    public:
      /** @name Constructors */
      //@{

      /// constructor from Clusters
      ClusteringResult(const std::string& name, const std::vector<Cluster>& clusters, double singleClusterSumOfSquares);

      //@}
      /** @name Getters */
      //@{

      /// returns the name
      std::string name() const;

      /// returns the total sum of squares
      double sumOfSquares() const;

      /// returns the r squared metric
      double rSquared() const;

      /// returns the adjusted r squared metric
      double rSquaredAdjusted() const;

      /// returns the current Clusters
      std::vector<Cluster> clusters() const;
      
      /// returns the next best ClusteringResult
      boost::optional<ClusteringResult> nextClusteringResult() const;

      //@}

    private:

      REGISTER_LOGGER("utilities.ClusteringResult");

      std::string m_name;
      std::vector<Cluster> m_clusters;
      double m_singleClusterSumOfSquares;
      double m_sumOfSquares;
      double m_rSquared;
      double m_rSquaredAdjusted;
  };

  /** AgglomerativeClusterer computes ClusteringResults and can help choose the best one.
  **/
  class UTILITIES_API AgglomerativeClusterer
  {
  public:
    /** @name Constructors */
    //@{

    /// constructor from TimeSeries, throws if TimeSeries does not contain equal number of samples per day
    AgglomerativeClusterer(const TimeSeries& timeSeries);

    //@}

    /// solves for the ClusteringResults, does nothing if there is already a solution
    bool solve();

    /// clears the current ClusteringResults
    void clear ();

    /// returns the ClusteringResults, must call solve before calling this
    std::vector<ClusteringResult> clusteringResults() const;

    /// returns the special ClusteringResults
    std::vector<ClusteringResult> specialClusteringResults() const;

  private:

    REGISTER_LOGGER("utilities.AgglomerativeClusterer");

    std::shared_ptr<ClusterData> m_clusterData;
    std::vector<ClusteringResult> m_clusteringResults;
    std::vector<ClusteringResult> m_specialClusteringResults; // computed in ctor
    double m_singleClusterSumOfSquares;
  };


} // openstudio

#endif // UTILITIES_DATA_AGGLOMERATIVECLUSTERING_HPP

