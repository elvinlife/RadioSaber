#include <cstdio>
#include <utility>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <vector>
using std::vector;

using coord_t = std::pair<int, int>;
using coord_cqi_t = std::pair<coord_t, int>;

static vector<int> MaximizeCell(int** flow_cqi, vector<int>& slice_quota_rbgs, int nb_rbgs, int nb_slices)
{
  // it's ok that the quota is negative
  vector<coord_cqi_t> sorted_cqi;
  vector<int>         slice_rbgs(nb_slices, 0);
  vector<int>         rbg_to_slice(nb_rbgs, -1);
  for (int i = 0; i < nb_rbgs; ++i)
    for (int j = 0; j < nb_slices; ++j) {
      sorted_cqi.emplace_back(coord_t(i, j), flow_cqi[i][j]);
    }
  std::sort(sorted_cqi.begin(), sorted_cqi.end(), [](coord_cqi_t a, coord_cqi_t b){return a.second > b.second; });
  for (auto it = sorted_cqi.begin(); it != sorted_cqi.end(); ++it) {
    int rbg_id = it->first.first;
    int slice_id = it->first.second;
    if (slice_rbgs[slice_id] < slice_quota_rbgs[slice_id] && rbg_to_slice[rbg_id] == -1) {
      rbg_to_slice[rbg_id] = slice_id;
      slice_rbgs[slice_id] += 1;
    } 
  }
  return rbg_to_slice;
}

static vector<int> VogelApproximate(int** flow_cqi, vector<int>& slice_quota_rbgs, int nb_rbgs, int nb_slices)
{
  vector<int>   slice_rbgs(nb_slices, 0);
  vector<int>   rbg_to_slice(nb_rbgs, -1);
  for (int i = 0; i < nb_rbgs; ++i) {
    int max_diff = -1;
    coord_t coord_1st;
    coord_t coord_2nd;
    // horizonal search
    for (int j = 0; j < nb_rbgs; ++j) {
      // the rbg has been allocated
      if (rbg_to_slice[j] != -1)
        continue;
      int cqi_1st = -1, cqi_2nd = -1;
      int slice_1st, slice_2nd;

      for (int k = 0; k < nb_slices; ++k) {
        // the slice has reached quota
        if (slice_rbgs[k] >= slice_quota_rbgs[k])
          continue;
        if (cqi_1st == -1 || flow_cqi[j][k] > cqi_1st) {
          slice_1st = k;
          cqi_1st = flow_cqi[j][k];
          continue;
        }
        if (cqi_2nd == -1 || flow_cqi[j][k] > cqi_2nd) {
          slice_2nd = k;
          cqi_2nd = flow_cqi[j][k];
          continue;
        }
      }
      if (cqi_1st - cqi_2nd > max_diff) {
        max_diff = cqi_1st - cqi_2nd;
        coord_1st = coord_t(j, slice_1st);
        coord_2nd = coord_t(j, slice_2nd);
      }
    }
    // vertical search
    for (int k = 0; k < nb_slices; ++k) {
      // the slice has reached quota
      if (slice_rbgs[k] >= slice_quota_rbgs[k])
        continue;
      int cqi_1st = -1, cqi_2nd = -1;
      int rbg_1st, rbg_2nd;
      for (int j = 0; j < nb_rbgs; ++j) {
        if (rbg_to_slice[j] != -1)
          continue;
        if (cqi_1st == -1 || flow_cqi[j][k] > cqi_1st) {
          rbg_1st = j;
          cqi_1st = flow_cqi[j][k];
          continue;
        }
        if (cqi_2nd == -1 || flow_cqi[j][k] > cqi_2nd) {
          rbg_2nd = j;
          cqi_2nd = flow_cqi[j][k];
          continue;
        }
      }
      if (cqi_1st - cqi_2nd > max_diff) {
        max_diff = cqi_1st - cqi_2nd;
        coord_1st = coord_t(rbg_1st, k);
        coord_2nd = coord_t(rbg_2nd, k);
      }
    }
    rbg_to_slice[coord_1st.first] = coord_1st.second;
    slice_rbgs[coord_1st.second] += 1;
  }
  return rbg_to_slice;
}

int main()
{
    int nb_rbgs = 5;
    int nb_slices = 3;
    int ** flow_cqi = new int*[nb_rbgs];
    for (int i = 0; i < nb_rbgs; ++i)
        flow_cqi[i] = new int[nb_slices];
    flow_cqi[0][0] = 8;
    flow_cqi[1][0] = flow_cqi[2][0] = flow_cqi[3][0] = flow_cqi[4][0] = 2;
    flow_cqi[0][1] = flow_cqi[1][1] = flow_cqi[2][1] = flow_cqi[3][1] = 10;
    flow_cqi[4][1] = 10;
    flow_cqi[0][2] = 9;
    flow_cqi[1][2] = 7;
    flow_cqi[2][2] = flow_cqi[3][2] = flow_cqi[4][2] = 2;
    vector<int> slice_quota_rbgs;
    slice_quota_rbgs.push_back(1);
    slice_quota_rbgs.push_back(3);
    slice_quota_rbgs.push_back(1);

    //auto rbg_to_slice = MaximizeCell(flow_cqi, slice_quota_rbgs, nb_rbgs, nb_slices); 
    auto rbg_to_slice = VogelApproximate(flow_cqi, slice_quota_rbgs, nb_rbgs, nb_slices); 
    for (auto it = rbg_to_slice.begin(); it != rbg_to_slice.end(); ++it)
        std::cout << *it << "; ";
    std::cout << std::endl;
}