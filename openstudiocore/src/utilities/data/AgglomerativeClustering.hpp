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

  class TimeSeries;

  /** A Cluster is a vector of Vectors that have been combined together.
  **/
  class UTILITIES_API Cluster
  {
  public:
    /** @name Constructors */
    //@{

    /// constructor from Vectors, throws if all vectors are not the same length or if length is 0
    Cluster(const std::vector<Vector>& vectors);

    //@}
    /** @name Getters */
    //@{

    /// returns the Vectors in this cluster
    std::vector<Vector> vectors() const;

    //@}

    /// computes the mean vector, $\hat x_{j} = \frac{1}{n} \sum_{i=0}^n x_{j}
    Vector meanVector() const;

    /// computes the within cluster sum of square error, $ss = \frac{1}{n} \sum_{i=0}^n \sum_{j=0}^m (x_{j}-\hat x_{j})^2
    double sumOfSquares() const;

    /// merges this cluster with another one, throws if all vectors are not the same length
    Cluster merge(const Cluster& other) const;

  private:

    REGISTER_LOGGER("utilities.Cluster");
    
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
      ClusteringResult(const std::string& name, double score, const std::vector<Cluster>& clusters);

      //@}
      /** @name Getters */
      //@{

      /// returns the name
      std::string name() const;

      /// returns the score
      double score() const;

      /// returns the Clusters
      std::vector<Cluster> clusters() const;

      //@}

    private:

      REGISTER_LOGGER("utilities.ClusteringResult");

      std::string m_name;
      double m_score;
      std::vector<Cluster> m_clusters;
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

    /// returns the ClusteringResults
    std::vector<ClusteringResult> clusterResults() const;

  private:

    REGISTER_LOGGER("utilities.AgglomerativeClusterer");

    std::vector<Vector> m_vectors;
    std::vector<ClusteringResult> m_clusterResults;
    std::vector<ClusteringResult> m_specialClusterResults; // computed in ctor
  };


} // openstudio

#endif // UTILITIES_DATA_AGGLOMERATIVECLUSTERING_HPP

