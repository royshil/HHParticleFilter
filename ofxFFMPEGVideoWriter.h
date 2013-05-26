//
//  ofxFFMPEGVideoWriter.h
//  HHParticleFilter
//
//  Created by roy_shilkrot on 4/7/13.
//
//  Copyright (c) 2013 MIT
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//  The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#ifndef __ofxFFMPEGVideoWriter__
#define __ofxFFMPEGVideoWriter__

#include <iostream>
#ifdef __cplusplus
extern "C" {
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
}
#endif

class ofxFFMPEGVideoWriter {
    //instance variables
    AVCodec *codec;
    int size, frame_count;
    AVFrame *picture, *picture_rgb24;
    struct SwsContext *sws_ctx;
    AVOutputFormat *fmt;
    AVFormatContext *oc;
    AVStream *video_st;
    AVCodecContext* c;

    bool initialized;

public:
    ofxFFMPEGVideoWriter():oc(NULL),codec(NULL),initialized(false),frame_count(1) {}
    
    /**
     * setup the video writer
     * @param output filename, the codec and format will be determined by it. (e.g. "xxx.mpg" will create an MPEG1 file
     * @param width of the frame
     * @param height of the frame
     * @param the bitrate
     * @param the framerate
     **/
    void setup(const char* filename, int width, int height, int bitrate = 400000, int framerate = 25);
    /**
     * add a frame to the video file
     * @param the pixels packed in RGB (24-bit RGBRGBRGB...)
     **/
    void addFrame(const uint8_t* pixels);
    /**
     * close the video file and release all datastructs
     **/
    void close();
    /**
     * is the videowriter initialized?
     **/
    bool isInitialized() const { return initialized; }
};

#endif /* defined(__ofxFFMPEGVideoWriter__) */
