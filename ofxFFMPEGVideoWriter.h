//
//  ofxFFMPEGVideoWriter.h
//  ShapeDeform
//
//  Created by roy_shilkrot on 4/7/13.
//
//

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
