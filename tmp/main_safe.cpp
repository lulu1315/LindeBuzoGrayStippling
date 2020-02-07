/*
 *      This is an interactive demo application for the algorithm proposed in:
 *
 *      Weighted Linde-Buzo Gray Stippling
 *      Oliver Deussen, Marc Spicker, Qian Zheng
 *
 *      In: ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia 2017)
 *      https://doi.org/10.1145/3130800.3130819
 *
 *     Copyright 2017 Marc Spicker (marc.spicker@googlemail.com)
 */
#include <cstdio>
#include <iostream>
#include <fstream>

#include "opencv2/highgui.hpp"
#include "opencv2/video.hpp"
#include "opencv2/core/ocl.hpp"
#include "opencv2/optflow.hpp"

#include <QApplication>
#include "mainwindow.h"
#include "lbgstippling.h"
#include <stdlib.h>

using namespace cv;
using namespace optflow;

int main(int argc, char* argv[]) {
    
    char outputimage[200];
    char outputply[200];
    std::ofstream outply;
    std::ifstream inply;
    float finalhysteresis;
    Mat_<Point2f> flow;
    
    LBGStippling::Params params;
    //stipples structure
    std::vector<Stipple> mystipples;
    std::vector<Stipple> initialstipples;
    //input args
    char *inputfile =                argv[1];
    char *outputroot =               argv[2];
    params.initialPoints=       atoi(argv[3]);
    params.initialPointSize=    atof(argv[4]);
    params.adaptivePointSize=   atoi(argv[5]);
    params.pointSizeMin=        atof(argv[6]);
    params.pointSizeMax=        atof(argv[7]);
    params.superSamplingFactor= atoi(argv[8]);
    params.maxIterations=       atoi(argv[9]);
    params.hysteresis=          atof(argv[10]);
    params.hysteresisDelta=     atof(argv[11]);
    float stipplesizefactor=    atof(argv[12]);
    int flowmode=               atoi(argv[13]);
    int hysteresis_strategy=    atoi(argv[14]);
    int curframe =              atoi(argv[15]);
    char *initialply =               argv[16];
    char *flowfile1 =                argv[17]; //.flo file or previous frame
    char *flowfile2 =                argv[18]; // current frame
    
    
    QApplication app(argc, argv); //collect args before otherwise atof doesn't work???
    //https://stackoverflow.com/questions/35566459/segfault-when-accessing-qapplicationarguments
    //MainWindow window;
    
    if (argc == 1)
        {
        printf("- gui mode\n");
        //QApplication app(argc, argv);
        MainWindow window;
        window.show();   
        return app.exec();
        }
    else
        {
        QImage density = QImage(inputfile);
        //print params
        printf("\n- cli mode :\n");;
        printf("--- inputimage.size %d %d\n",density.width(),density.height());
        printf("--- params.initialPoints %zu\n",params.initialPoints);
        printf("--- params.initialPointSize %f\n",params.initialPointSize);
        printf("--- params.adaptivePointSize %d\n",params.adaptivePointSize);
        printf("--- params.pointSizeMin %f\n",params.pointSizeMin);
        printf("--- params.pointSizeMax %f\n",params.pointSizeMax);
        printf("--- params.superSamplingFactor %zu\n",params.superSamplingFactor);
        printf("--- params.maxIterations %zu\n",params.maxIterations);
        printf("--- params.hysteresis %f\n",params.hysteresis);
        printf("--- params.hysteresisDelta %f\n",params.hysteresisDelta);
        printf("--- stipplesizefactor %f\n",stipplesizefactor);
        printf("--- optical flow mode %d\n",flowmode);
        printf("--- hysteresis_strategy %d\n\n",hysteresis_strategy);
        printf("--- current frame %d\n\n",curframe);
        if (flowmode == 0) //no opticalflow
            {
            printf ("- no optical flow\n");
            //performs stippling
            LBGStippling m_stippling;
            mystipples = m_stippling.stipple(density,flow,params,0,flowmode,initialstipples,finalhysteresis,hysteresis_strategy,curframe);
            }
        if (flowmode >= 1) 
        //1: backward flow 2: forward flow 
        // assumes flowfile1 is a precomputed optical flow file in .flo format
        // backward_current_prev.flo or forward_prev_current.flo
        //3: backward flow 4: forward flow 
        // assumes flowfile1 is the previous pixel image
        // assumes flowfile2 is the current pixel image
        // calculate optical flow with deepflow
            {
            //read initial ply
            std::string plyword;
            printf ("- reading previous stipples result\n");
            inply.open (argv[16]);
            //find previous finalhysteresis
            while (plyword != "finalhysteresis")
                {
                inply >> plyword;
                }
            inply >> finalhysteresis;
            printf("--- previous finalhysteresis %f\n",finalhysteresis);
            //find number of vertices
            while (plyword != "vertex")
                {
                inply >> plyword;
                }
            int nvertices;
            inply >> nvertices;
            //printf("%d\n",nvertices);
            //stipples structures
            std::vector<Stipple> initialstipples;
            //continue in ply
            while (plyword != "end_header")
                {
                inply >> plyword;
                }
            //fill initialstipples with ply values
            float init_stipple_xpos,init_stipple_ypos,init_stipple_zpos,init_stipple_xflow,init_stipple_zflow,init_stipple_size;
            int init_stipple_label,init_stipple_birth;
            for (int i=1;i<=nvertices;i++)
                {
                inply >> init_stipple_xpos;
                inply >> init_stipple_ypos;
                inply >> init_stipple_zpos;
                inply >> init_stipple_xflow;
                inply >> init_stipple_zflow;
                inply >> init_stipple_size;
                inply >> init_stipple_label;
                inply >> init_stipple_birth;
                //printf("%f %f %f %f\n",init_stipple_xpos,init_stipple_ypos,init_stipple_zpos,init_stipple_size);
                initialstipples.push_back({QVector2D(init_stipple_xpos,init_stipple_zpos),QVector2D(init_stipple_xflow,init_stipple_zflow), init_stipple_size, init_stipple_label, init_stipple_birth, Qt::black});
                }
            inply.close();
            printf("--- initial stipples size : %zu\n",initialstipples.size());
            
            //read precomputed flow file
            if (flowmode == 1 || flowmode == 2)
                {
                if (flowmode == 1) printf ("\n- using backward precomputed optical flow : %s\n",flowfile1);
                if (flowmode == 2) printf ("\n- using forward precomputed optical flow : %s\n",flowfile1);
                flow = optflow::readOpticalFlow(flowfile1);
                printf("--- flow size : %d %d\n",flow.size[1], flow.size[0]);
                }
            //compute optical flow on the fly
            if (flowmode == 3 || flowmode == 4)
                {
                // read current and previous frames in opencv format
                //should be better to convert density to cv::Mat ...
                Mat cur, prev;
                prev= imread(flowfile1, 1);
                cur = imread(flowfile2, 1);
                //convert to grayscale
                if ( cur.depth() != CV_8U ) cur.convertTo(cur, CV_8U);
                if ( prev.depth() != CV_8U ) prev.convertTo(prev, CV_8U);
                cvtColor(cur, cur, COLOR_BGR2GRAY);
                cvtColor(prev, prev, COLOR_BGR2GRAY);
                //compute opencv deep optical flow
                flow = Mat(cur.size[0], cur.size[1], CV_32FC2);
                Ptr<DenseOpticalFlow> algorithm;
                algorithm = createOptFlow_DeepFlow();
                //use openCL
                int useGpu = 1;
                cv::ocl::setUseOpenCL(useGpu);
                //go into the flow
                if (flowmode == 3) //backward
                    {
                    printf ("\n- computing backward optical flow\n");
                    printf ("--- OpenCL Enabled: %u\n", useGpu && cv::ocl::haveOpenCL());
                    if (useGpu) algorithm->calc(cur, prev, flow.getUMat(ACCESS_RW));
                    else algorithm->calc(cur, prev, flow);
                    }
                if (flowmode == 4) //forward
                    {
                    printf ("\n- computing forward optical flow\n");
                    printf ("--- OpenCL Enabled: %u\n", useGpu && cv::ocl::haveOpenCL());
                    if (useGpu) algorithm->calc(prev, cur, flow.getUMat(ACCESS_RW));
                    else algorithm->calc(prev, cur, flow);
                    }
                }
            //warp initialstipples with OpticalFlow
            for (auto &stipple : initialstipples) {
                float pixelflowx,pixelflowy;
                pixelflowx=stipple.pos[0]*density.width();
                pixelflowy=stipple.pos[1]*density.height();
                Vec2f of = flow.at<Vec2f>(pixelflowy,pixelflowx);
                //printf("flow (x,y) : (%f,%f) %f %f\n",pixelflowx,pixelflowy,of.val[0],of.val[1]);
                if (flowmode == 1 || flowmode == 3) {
                stipple.pos[0] = stipple.pos[0] - (of.val[0]/(float)density.width());
                stipple.pos[1] = stipple.pos[1] - (of.val[1]/(float)density.height());}
                if (flowmode == 2 || flowmode == 4) {
                stipple.pos[0] = stipple.pos[0] + (of.val[0]/(float)density.width());
                stipple.pos[1] = stipple.pos[1] + (of.val[1]/(float)density.height());}
                }
            //performs stippling
            LBGStippling m_stippling;
            mystipples = m_stippling.stipple(density,flow,params,2,flowmode,initialstipples,finalhysteresis,hysteresis_strategy,curframe);
            }
        //prepare drawings
        QImage background=density.copy();
        background.fill(qRgb(255, 255, 255));
        QPainter painter(&background);
        QBrush   brush(Qt::black);
        QPen     pen(Qt::NoPen);
        painter.setPen(pen);
        painter.setBrush(brush);
        //printf("--- finalhysteresis %f\n",finalhysteresis);
        //ply 
        sprintf(outputply,"%s.ply",outputroot);
        outply.open (outputply);
        outply << "ply\n";
        outply << "format ascii 1.0\n";
        outply << "comment finalhysteresis " << finalhysteresis << "\n";
        outply << "comment created by LBGStippling\n";
        outply << "element vertex " << mystipples.size() << "\n";
        outply << "property float x\n";
        outply << "property float y\n";
        outply << "property float z\n";
        outply << "property float xflow\n";
        outply << "property float zflow\n";
        outply << "property float size\n";
        outply << "property int label\n";
        outply << "property int birth\n";
        outply << "element face 0\n";
        outply << "property list uchar int vertex_indices\n";
        outply << "end_header\n";
        //draw stipples anf fill .ply file
        for (const auto &stipple : mystipples) {
            //colored stipples from image
            //QRgb pixcolor=density.pixel(stipple.pos[0]*density.width(),stipple.pos[1]*density.height());
            //brush.setColor(QColor(qRed(pixcolor),qGreen(pixcolor),qBlue(pixcolor)));
            //pen.setColor(QColor(qRed(pixcolor),qGreen(pixcolor),qBlue(pixcolor)));
            brush.setColor(stipple.color);
            //pen.setColor(stipple.color);
            painter.setPen(pen);
            painter.setBrush(brush);
            painter.setRenderHint(QPainter::Antialiasing, true);
            painter.drawEllipse(QPointF(stipple.pos[0]*density.width(),stipple.pos[1]*density.height()), stipple.size/stipplesizefactor, stipple.size/stipplesizefactor);
            outply << stipple.pos[0] << " 0 " << stipple.pos[1] << " " << stipple.flow[0] << " " << stipple.flow[1] << " " << stipple.size << " " << stipple.label << " " << stipple.birth << "\n";
            }
        painter.end();
        outply.close();
        //save result image
        sprintf(outputimage,"%s.png",outputroot);
        printf ("\n- writing image      : %s\n",outputimage);
        printf ("- writing pointcloud : %s\n",outputply);
        background.save(outputimage);
        return 1;
        }
}
