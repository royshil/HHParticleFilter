/*
 *  main,cpp
 *  HHParticleFilter
 *
 *  Created by Roy Shilkrot on 03/15/13.
 *
 *  Example usage of the Heap & Hogg Particle Filter for hand curve tracking.
 *
 *  Copyright (c) 2013 MIT
 *
 *  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *  The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 */


#include "std.h"

using namespace cv;

#include "CurveCSS.h"
#include "CurveSignature.h"

#include "HPDMTracker.h"
#include "ofxFFMPEGVideoWriter.h"


Mat cameraframe;
HPDMTracker tracker;


cv::Point start_touch,last_touch;
void onMouse(int e, int x, int y, int flags, void* ) {
	static Point2d start_touch(-1,-1);
	static Point2d last_touch(-1,-1);
    	
	Mat show; cameraframe.copyTo(show);
    putText(show, "click on center and drag to scale", cv::Point(10,show.rows-20), CV_FONT_HERSHEY_PLAIN, 1.0, Scalar(0,0,255));
    
	Point2d touch(x,y);
	Point2d d = start_touch-touch;
    if( e == CV_EVENT_LBUTTONDOWN ) {
		//press
		start_touch = last_touch = touch;
	} else if( e == CV_EVENT_LBUTTONUP ) {
		//release
		tracker.setInitialGuess(cv::Point(start_touch.x,start_touch.y), norm(d)/200.0, atan2(d.y, -d.x) / CV_PI * 180.0);
        tracker.visualizeParticles(show);
	} else if (e == CV_EVENT_MOUSEMOVE && flags == CV_EVENT_FLAG_LBUTTON) {
		//drag
		line(show, start_touch, touch, Scalar(0,0,255), 1);
		stringstream ss; ss << "center " << start_touch << " scale " << norm(d)/200.0 << " angle " << (atan2(d.y, -d.x) / CV_PI * 180.0);
		putText(show, ss.str(), cv::Point(10,15), CV_FONT_HERSHEY_PLAIN, 1.0, Scalar(0,0,255), 1);
		
		last_touch = touch;
	}
	
	imshow("tracker", show);
	waitKey(30);
}


int main(int argc, char** argv) {
    VideoCapture capture("recording_mask.avi");
	if (!capture.isOpened()) {
		cerr << "cannot open video\n";
		exit(0);
	}
    
    //    tracker.initializeFromDirectory("training_curves");
    //	destroyAllWindows();
    //    tracker.saveToFile("HPDMTracker.xml");
    tracker.initializeFromFile("HPDMTracker.xml",true);
    destroyAllWindows();

    capture >> cameraframe;
    namedWindow("tracker");
	setMouseCallback("tracker", onMouse, 0);
    vector<Point2f> mycurve;
    GetCurveForImage(cameraframe, mycurve, false);
    SimpleSmoothCurve(mycurve, mycurve, 3.0);
    drawOpenCurve(cameraframe, mycurve, Scalar(255), 2);
	imshow("tracker", cameraframe);
	waitKey(0);
    
    ofxFFMPEGVideoWriter writer;

    Mat tmp,tmp_small;
	for (;;) {
        capture >> cameraframe;
        if (cameraframe.empty()) {
            break;
        }
		
		GaussianBlur(cameraframe, cameraframe, Size(11,11), 3.0, 3.0);
		threshold(cameraframe, cameraframe, 128.0, 255.0, CV_THRESH_BINARY);
		
		tracker.update(cameraframe,Mat());

		Mat show; cameraframe.copyTo(show);
		drawOpenCurve(show, tracker.getMeanShape(), Scalar(0,0,255), 1);
		tracker.visualizeParticles(show);
        
		imshow("tracker", show);
        if(!writer.isInitialized())
            writer.setup("output.gif", show.cols/2, show.rows/2);
        {
            cvtColor(show, tmp, CV_BGR2RGB);
            resize(tmp, tmp_small, Size(show.cols/2,show.rows/2));
            writer.addFrame(tmp_small.data);
        }
        int c = -1;
		if((c = waitKey(15))==27) break;
        if(c == ' ') waitKey();
	}
    writer.close();
}
