#include <ceres/rotation.h>

template<class T>
bool ProjectionErrorFunction::operator()(const T* const camera,
                                         const T* const object_point,
                                         T* error) const {
  // Rotate into camera co-ordinates.
  T camera_point[3];
  ceres::UnitQuaternionRotatePoint(camera, object_point, camera_point);

  // Subtract from point.
  for (int d = 0; d < 2; d += 1) {
    error[d] = w_[d] - camera_point[d];
  }

  return true;
}
