#include "lbgstippling.h"
#include "voronoicell.h"

#include <cassert>
#include <random>

#include <QVector>
#include <QtMath>

#include <cstdio>
#include <iostream>
#include <fstream>

namespace Random {
//make it reproductible
//static std::random_device rd;
//static std::mt19937 gen(rd());
//https://diego.assencio.com/?index=6890b8c50169ef45b74db135063c227c
float seed=33.33;
static std::mt19937 gen(seed);
}  // namespace Random

using Params = LBGStippling::Params;
using Status = LBGStippling::Status;

QVector<QVector2D> sites(const std::vector<Stipple> &stipples) {
  QVector<QVector2D> sites(stipples.size());
  std::transform(stipples.begin(), stipples.end(), sites.begin(),
                 [](const auto &s) { return s.pos; });
  return sites;
}

std::vector<Stipple> randomStipples(size_t n, float size) {
  std::uniform_real_distribution<float> dis(0.01f, 0.99f);
  std::vector<Stipple> stipples(n);
  std::generate(stipples.begin(), stipples.end(), [&]() {
    return Stipple{QVector2D(dis(Random::gen),dis(Random::gen)), QVector2D(0,0), QVector2D(0,0) , size,
                   Qt::black};
  });
  return stipples;
}
    
template <class T>
inline T pow2(T x) {
  return x * x;
}

QVector2D jitter(QVector2D s) {
  using namespace Random;
  std::uniform_real_distribution<float> jitter_dis(-0.001f, 0.001f);
  return s += QVector2D(jitter_dis(gen), jitter_dis(gen));
}

float getSplitValueUpper(float pointDiameter, float hysteresis,
                         size_t superSampling) {
  const float pointArea = M_PIf32 * pow2(pointDiameter / 2.0f);
  return (1.0f + hysteresis / 2.0f) * pointArea * pow2(superSampling);
}

float getSplitValueLower(float pointDiameter, float hysteresis,
                         size_t superSampling) {
  const float pointArea = M_PIf32 * pow2(pointDiameter / 2.0f);
  return (1.0f - hysteresis / 2.0f) * pointArea * pow2(superSampling);
}

float stippleSize(const VoronoiCell &cell, const Params &params) {
  if (params.adaptivePointSize) {
    const float avgIntensitySqrt = std::sqrt(cell.sumDensity / cell.area);
    return params.pointSizeMin * (1.0f - avgIntensitySqrt) +
           params.pointSizeMax * avgIntensitySqrt;
  } else {
    return params.initialPointSize;
  }
}

float currentHysteresis(size_t i, const Params &params) {
  return params.hysteresis + i * params.hysteresisDelta;
}

float updateHysteresis(size_t i, float hysteresis,float hysteresisDelta) {
  return hysteresis + i * hysteresisDelta;
}

bool notFinished(const Status &status, const Params &params) {
  auto [iteration, size, splits, merges, hysteresis] = status;
  return !((splits == 0 && merges == 0) || (iteration == params.maxIterations));
}

LBGStippling::LBGStippling() {
  m_statusCallback = [](const Status &) {};
  m_stippleCallback = [](const std::vector<Stipple> &) {};
}

void LBGStippling::setStatusCallback(Report<Status> statusCB) {
  m_statusCallback = statusCB;
}

void LBGStippling::setStippleCallback(Report<std::vector<Stipple>> stippleCB) {
  m_stippleCallback = stippleCB;
}

std::vector<Stipple> LBGStippling::stipple(const QImage &density,
                                           const Mat_<Point2f> &previousflow,
                                           const Mat_<Point2f> &nextflow,
                                           const Params &params,
                                           const int &mode,
                                           const int &flowmode,
                                           const std::vector<Stipple> &initialstipples,
                                           float &finalhysteresis,
                                           const int &hysteresis_strategy,
                                           const int &curframe) const {
  QImage densityGray =
      density
          .scaledToWidth(params.superSamplingFactor * density.width(),
                         Qt::SmoothTransformation)
          .convertToFormat(QImage::Format_Grayscale8);
          
  VoronoiDiagram voronoi(densityGray);
  std::vector<Stipple> stipples;
  
  if (mode == 0 || mode == 1)
    {
    stipples = randomStipples(params.initialPoints, params.initialPointSize);
    }
    
  if (mode == 2 || mode == 3)
    {
    stipples = initialstipples;
    }
    
  Status status = {0, 0, 1, 1, params.hysteresis};
  
  printf("\n- stippling : iterations : ");
  while (notFinished(status, params)) {
    
    status.splits = 0;
    status.merges = 0;
    auto indexMap = voronoi.calculate(sites(stipples));
    std::vector<VoronoiCell> cells = accumulateCells(indexMap, densityGray);

    //transfer labels and birth to cells
    int maxlabel=0;
    for (int k = 0; k < cells.size(); ++k) {
        cells[k].label=stipples[k].label;
        cells[k].birth=stipples[k].birth;
        if (stipples[k].label > maxlabel) {maxlabel=stipples[k].label;}
        }
        
    /*
    //looking at cells ...
    QImage densityCheck(indexMap.width,indexMap.height,QImage::Format_RGB32);
    for (int x = 0; x < indexMap.width; ++x) {
        for (int y = 0; y < indexMap.height; ++y) {
        uint32_t index = indexMap.get(x, y);
        VoronoiCell& checkcell = cells[index];
        float avgIntensitySqrt = std::sqrt(checkcell.sumDensity / checkcell.area);
        if (avgIntensitySqrt > 1.) {printf("x/y >1 : %d/%d : %f\n",x,y,avgIntensitySqrt);}
        //float avgIntensitySqrt = checkcell.sumDensity / checkcell.area;
        densityCheck.setPixelColor(x,y,QColor(avgIntensitySqrt*255.0,avgIntensitySqrt*255.0,avgIntensitySqrt*255.0));
        }
    }
    char tmpout[200];
    sprintf(tmpout,"tmp.%zu.png",status.iteration+1);
    printf ("\ncheckimage : %s\n",tmpout);
    densityCheck.save(tmpout);
    //end looking at cells
    */
    
    //checking if cells are ordered like stipples...
    //for (int k = 0; k < stipples.size(); ++k) {
    //    printf("k : %d cell centroid : %f %f stipple pos : %f %f\n",k,cells[k].centroid.x(),cells[k].centroid.y(),stipples[k].pos.x(),stipples[k].pos.y());
     //   }
      
    assert(cells.size() == stipples.size());

    stipples.clear();
    float hysteresis;
    
    if (mode == 0)
    {
    hysteresis = currentHysteresis(status.iteration, params);
    }
    if (mode == 1 || mode == 2 || mode == 3 || mode == 4)
    {
        if (hysteresis_strategy == 0)
            {
            hysteresis = currentHysteresis(status.iteration, params);
            }
        if (hysteresis_strategy == 1)
            {
            //hysteresis = finalhysteresis;
            hysteresis = updateHysteresis(status.iteration,finalhysteresis,params.hysteresisDelta);
            }
        if (hysteresis_strategy >= 2)
            {
            hysteresis = updateHysteresis(status.iteration, (((hysteresis_strategy-1)*finalhysteresis)+params.hysteresis)/(hysteresis_strategy), params.hysteresisDelta);
            }
    }
    
    status.hysteresis = hysteresis;
    int count = 0;
    for (const auto &cell : cells) {
      const float totalDensity = cell.sumDensity;
      float diameter = stippleSize(cell, params);
      
      ++count;
      
      if (totalDensity < getSplitValueLower(diameter, hysteresis,
                                            params.superSamplingFactor) ||
          cell.area == 0.0f) {
        // cell too small - merge
        ++status.merges;
        //printf("-%d",cell.label);
        continue;
      }

      if (totalDensity < getSplitValueUpper(diameter, hysteresis,
                                            params.superSamplingFactor)) {
        // cell size within acceptable range - keep
        stipples.push_back({cell.centroid, QVector2D(0,0), QVector2D(0,0) ,0, diameter, cell.label, cell.birth, Qt::black});
        continue;
      }

      //printf("index : %d totaldensity : %f Lu/Ld %f %f\n",count,totalDensity,getSplitValueUpper(diameter, hysteresis,1),getSplitValueLower(diameter, hysteresis,1));
      // cell too large - split
      const float area = std::max(1.0f, cell.area);
      const float circleRadius = std::sqrt(area / M_PIf32);
      QVector2D splitVector = QVector2D(0.5f * circleRadius, 0.0f);

      const float a = cell.orientation;
      QVector2D splitVectorRotated = QVector2D(
          splitVector.x() * std::cos(a) - splitVector.y() * std::sin(a),
          splitVector.y() * std::cos(a) + splitVector.x() * std::sin(a));

      splitVectorRotated.setX(splitVectorRotated.x() / densityGray.width());
      splitVectorRotated.setY(splitVectorRotated.y() / densityGray.height());

      QVector2D splitSeed1 = cell.centroid - splitVectorRotated;
      QVector2D splitSeed2 = cell.centroid + splitVectorRotated;

      // check boundaries
      splitSeed1.setX(std::max(0.0f, std::min(splitSeed1.x(), 1.0f)));
      splitSeed1.setY(std::max(0.0f, std::min(splitSeed1.y(), 1.0f)));

      splitSeed2.setX(std::max(0.0f, std::min(splitSeed2.x(), 1.0f)));
      splitSeed2.setY(std::max(0.0f, std::min(splitSeed2.y(), 1.0f)));

      //stipples.push_back({jitter(splitSeed1), diameter, Qt::red});
      //stipples.push_back({jitter(splitSeed2), diameter, Qt::red});
      //maxlabel ++;
      stipples.push_back({splitSeed1, QVector2D(0,0), QVector2D(0,0) ,0, diameter, cell.label , cell.birth , Qt::red});
      maxlabel ++;
      stipples.push_back({splitSeed2, QVector2D(0,0), QVector2D(0,0) ,0, diameter, maxlabel , curframe , Qt::red});
      //printf("+%d",maxlabel);

      ++status.splits;
    }
    
    status.size = stipples.size();
    m_stippleCallback(stipples);
    m_statusCallback(status);

    ++status.iteration;
    printf(".");
    fflush(stdout);
  }
  
  printf("\n--- done %zu iterations\n",status.iteration);
  printf("--- status.size %zu\n",status.size);
  printf("--- status.splits %zu\n",status.splits);
  printf("--- status.merges %zu\n",status.merges);
  printf("--- status.hysteresis %f\n",status.hysteresis);
  finalhysteresis = status.hysteresis;
    
  //motion vectors evaluation
  if (flowmode == 1) 
    {
    printf("--- evaluating motion vectors\n");
    //recreate indexmap with finalstipples
    auto indexMap = voronoi.calculate(sites(stipples));
    std::vector<Vec2f> averagepreviousof = std::vector<Vec2f>(indexMap.count());
    std::vector<Vec2f> averagenextof = std::vector<Vec2f>(indexMap.count());
    std::vector<float> luminances = std::vector<float>(indexMap.count());
    std::vector<int>   cellcount = std::vector<int>(indexMap.count());
    for (int x = 0; x < indexMap.width; ++x) {
        for (int y = 0; y < indexMap.height; ++y) {
            uint32_t index = indexMap.get(x, y);
            QRgb densityPixel = density.pixel(x, y);
            float luminance = qGray(densityPixel);
            Vec2f prevof = previousflow.at<Vec2f>(y,x);
            Vec2f nextof = nextflow.at<Vec2f>(y,x);
            Vec2f& accumulateprevof = averagepreviousof[index];
            Vec2f& accumulatenextof = averagenextof[index];
            float& accumulateluminance = luminances[index];
            int& divider            = cellcount[index];
            divider++;
            accumulateprevof=accumulateprevof+prevof;
            accumulatenextof=accumulatenextof+nextof;
            accumulateluminance=accumulateluminance+luminance;
            }
        }
    for (int k = 0; k < stipples.size(); ++k) {
        //printf("cell : %d cellcount : %d of : %f/%f\n",k,cellcount[k],averageof[k][0]/cellcount[k],averageof[k][1]/cellcount[k]);
        stipples[k].previousflow=QVector2D(averagepreviousof[k][0]/cellcount[k],averagepreviousof[k][1]/cellcount[k]);
        stipples[k].nextflow=QVector2D(averagenextof[k][0]/cellcount[k],averagenextof[k][1]/cellcount[k]);
        stipples[k].luminance=luminances[k]/cellcount[k];
        }
    }
    
  if (mode == 0 || mode == 1)
    {
    //create labels
    for (int k = 0; k < stipples.size(); ++k) {
        stipples[k].label=k;
        stipples[k].birth=curframe;
        }
    }
    
    int maxlabel=0;
    for (int k = 0; k < stipples.size(); ++k) {
        if (stipples[k].label > maxlabel) {maxlabel=stipples[k].label;}
        }
    printf("--- maxlabel %d\n",maxlabel);
    
  return stipples;
}
