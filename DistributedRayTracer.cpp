#include "ParallelRayTracer.h"
#include "imager.h"
void MPI_TraceRay(int& rank, int& numProc, Imager::PixelData* localBUffer, Imager::Scene* scene,
        size_t pixelsWide,
        size_t pixelsHigh,
        double zoom,
        size_t antiAliasFactor)
{
    using namespace Imager;
    // Oversample the image using the anti-aliasing factor.
    const size_t largePixelsWide = antiAliasFactor * pixelsWide;
    const size_t largePixelsHigh = antiAliasFactor * pixelsHigh;
    const size_t smallerDim =
        ((pixelsWide < pixelsHigh) ? pixelsWide : pixelsHigh);

    const double largeZoom = antiAliasFactor * zoom * smallerDim;
    // ImageBuffer buffer(largePixelsWide, largePixelsHigh, backgroundColor);
    // The camera is located at the origin.
    Vector camera(0.0, 0.0, 0.0);
    Vector direction(0.0, 0.0, -1.0);

    const Color fullIntensity(1.0, 1.0, 1.0);
    pthread_t threads[THREAD_COUNT];
    printf("starting pthread\n");
    rt_chunk* tasks = new rt_chunk[THREAD_COUNT];
    int dx, dy, xt;
    xt = 0;
    dx = largePixelsWide / THREAD_COUNT;
    for (size_t i = 0; i < THREAD_COUNT; i++)
    {
        tasks[i].xStart = xt;
        tasks[i].xEnd = xt += dx;
        tasks[i].yStart = rank*(largePixelsHigh/numProc);
        tasks[i].yEnd = (rank+1)*(largePixelsHigh/numProc);
        tasks[i].scene = scene;
        tasks[i].refractiveIndex = scene->GetAmbientRefraction();
        tasks[i].rayIntensity = fullIntensity;
        tasks[i].vantage = &camera;
        tasks[i].smallerDim = smallerDim;
        tasks[i].largeZoom = largeZoom;
        tasks[i].largePixelsHigh = largePixelsHigh;
        tasks[i].largePixelsWide = largePixelsWide;
        tasks[i].lbuffer = localBUffer;
        tasks[i].imagewidth = pixelsWide;
        tasks[i].offset = largePixelsWide*(largePixelsHigh/numProc)*rank;
    }
    for (size_t i = 0; i < THREAD_COUNT; i++)
    {
        pthread_create(&threads[i], NULL, mt_TraceRay, (void*)&tasks[i]);
    }
    printf("waiting barrier\n");
    for (size_t i = 0; i < THREAD_COUNT; i++)
    {
        pthread_join(threads[i], NULL);
    }
    printf("done\n");
    printf("proc %d done\n", rank);
}

