CXXFLAGS=-O3
all: qcluster.cc freqs.cc  freqs.hh dists.hh dists.cc centroid.hh centroid.cc em.cc  em.hh  seqio.cc  seqio.hh
	g++ ${CXXFLAGS} -g -Wall qcluster.cc em.cc dists.cc centroid.cc freqs.cc seqio.cc -o qCluster

qcluster: qcluster.cc freqs.cc  freqs.hh dists.hh dists.cc centroid.hh centroid.cc em.cc  em.hh  seqio.cc  seqio.hh
	g++ ${CXXFLAGS} -g -Wall qcluster.cc em.cc dists.cc centroid.cc freqs.cc seqio.cc -o qCluster
