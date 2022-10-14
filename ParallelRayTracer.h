#include <cmath>
#include <fstream>
#include <iostream>
// #include <pthread.h>
#include "imager.h"
#include "lodepng.h"
#define THREAD_COUNT 4
#define MT_MAX_OPTICAL_RECURSION_DEPTH 5
struct rt_chunk
{
    int xStart;
    int xEnd;
    int yStart;
    int yEnd;
    Imager::PixelData* lbuffer;
    Imager::ImageBuffer* buffer;
    Imager::Vector* vantage;
    double refractiveIndex;
    Imager::Color rayIntensity;
    const Imager::Scene* scene;
    int smallerDim;
    int largeZoom;
    int largePixelsWide;
    int largePixelsHigh;
    int imagewidth;
    int offset;
};

void *mt_TraceRay(void* args);

void MPI_TraceRay(int& rank, int& numProc, Imager::PixelData* localBUffer, Imager::Scene* scene,
        size_t pixelsWide,
        size_t pixelsHigh,
        double zoom,
        size_t antiAliasFactor);

void SaveImageFromBuffer(
    const char *outPngFileName, 
    size_t pixelsWide, 
    size_t pixelsHigh, 
    double zoom, 
    size_t antiAliasFactor,
    Imager::ImageBuffer* _buffer);

static unsigned char ConvertPixelValue(
    double colorComponent, 
    double maxColorValue)
{
    int pixelValue = 
        static_cast<int> (255.0 * colorComponent / maxColorValue);

    // Clamp to the allowed range of values 0..255.
    if (pixelValue < 0)
    {
        pixelValue = 0;
    }
    else if (pixelValue > 255)
    {
        pixelValue = 255;
    }

    return static_cast<unsigned char>(pixelValue);
}