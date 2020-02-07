#ifndef LBGSTIPPLING_H
#define LBGSTIPPLING_H

#include "opencv2/highgui.hpp"
#include "opencv2/video.hpp"
#include "opencv2/core/ocl.hpp"
#include "opencv2/optflow.hpp"

#include "voronoidiagram.h"

#include <QImage>
#include <QVector2D>

using namespace cv;
using namespace optflow;

// TODO: Color is only used for debugging
struct Stipple {
  QVector2D pos;
  QVector2D previousflow;
  QVector2D nextflow;
  float luminance;
  float size;
  int label;
  int birth;
  QColor color;
};

class LBGStippling {
 public:
  struct Params {
    size_t initialPoints = 1;
    float initialPointSize = 4.0f;

    bool adaptivePointSize = true;
    float pointSizeMin = 2.0f;
    float pointSizeMax = 4.0f;

    size_t superSamplingFactor = 1;
    size_t maxIterations = 50;

    float hysteresis = 0.6f;
    float hysteresisDelta = 0.01f;
  };

  struct Status {
    size_t iteration;
    size_t size;
    size_t splits;
    size_t merges;
    float hysteresis;
  };

  template <class T>
  using Report = std::function<void(const T&)>;

  LBGStippling();

  std::vector<Stipple> stipple(const QImage& density,
                               const Mat_<Point2f>& previousflow,
                               const Mat_<Point2f>& nextflow,
                               const Params& params,
                               const int& mode,
                               const int& flowmode,
                               const std::vector<Stipple>& initialstipples,
                               float& finalhysteresis,
                               const int& hysteresis_strategy,
                               const int& curframe) const;

  // TODO: Rename and method chaining.
  void setStatusCallback(Report<Status> statusCB);
  void setStippleCallback(Report<std::vector<Stipple>> stippleCB);

 private:
  Report<Status> m_statusCallback;
  Report<std::vector<Stipple>> m_stippleCallback;
};

#endif  // LBGSTIPPLING_H
