#include <iostream>
#include "block.h"
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include "chessboard.h"
#include </usr/include/mpi/mpi.h>
#include "ParallelRayTracer.h"

#pragma warning(disable : 4996)
#define WIDTH 1920
#define HEIGHT 1080
#define ZOOM 1
#define AA_FACTOR 3
void SaveImageFromBuffer(
    const char *outPngFileName, 
    size_t pixelsWide, 
    size_t pixelsHigh, 
    double zoom, 
    size_t antiAliasFactor,
    Imager::ImageBuffer* _buffer)
    {
        const double max = _buffer->MaxColorValue();

        // Downsample the image buffer to an integer array of RGBA 
        // values that LodePNG understands.
        const unsigned char OPAQUE_ALPHA_VALUE = 255;
        const unsigned BYTES_PER_PIXEL = 4;

        // The number of bytes in buffer to be passed to LodePNG.
        const unsigned RGBA_BUFFER_SIZE =
            pixelsWide * pixelsHigh * BYTES_PER_PIXEL;

        std::vector<unsigned char> rgbaBuffer(RGBA_BUFFER_SIZE);
        unsigned rgbaIndex = 0;
        const double patchSize = antiAliasFactor * antiAliasFactor;
        for (size_t j = 0; j < pixelsHigh; ++j)
        {
            for (size_t i = 0; i < pixelsWide; ++i)
            {
                Imager::Color sum(0.0, 0.0, 0.0);
                for (size_t di = 0; di < antiAliasFactor; ++di)
                {
                    for (size_t dj = 0; dj < antiAliasFactor; ++dj)
                    {
                        sum += _buffer->Pixel(
                            antiAliasFactor * i + di,
                            antiAliasFactor * j + dj).color;
                    }
                }
                sum /= patchSize;

                // Convert to integer red, green, blue, alpha values,
                // all of which must be in the range 0..255.
                rgbaBuffer[rgbaIndex++] = ConvertPixelValue(sum.red, max);
                rgbaBuffer[rgbaIndex++] = ConvertPixelValue(sum.green, max);
                rgbaBuffer[rgbaIndex++] = ConvertPixelValue(sum.blue, max);
                rgbaBuffer[rgbaIndex++] = OPAQUE_ALPHA_VALUE;
            }
        }

        // Write the PNG file
        const unsigned error = lodepng::encode(
            outPngFileName,
            rgbaBuffer,
            pixelsWide,
            pixelsHigh);

        // If there was an encoding error, throw an exception.
        if (error != 0)
        {
            std::string message = "PNG encoder error: ";
            message += lodepng_error_text(error);
            throw Imager::ImagerException(message.c_str());
        }
    }

void constructLocalScene(Imager::Scene* scene)
{
    using namespace Imager;
    ChessBoard* board = new ChessBoard(50,1,1,.25, Color(0.75, 0.70, 0.10),Color(0.30, 0.30, 0.40), Color(0.50, 0.30, 0.10));
    board->Move(-0.35, -4, -20.0);
    board->RotateZ(+10.0);
    board->RotateX(-90);
    scene->AddSolidObject(board);
    ChessBoard* board2 = new ChessBoard(150, 1,  1, .25, Color(0.75, 0.70, 0.10), Color(0.30, 0.30, 0.40), Color(1, 1, 1));
    board2->Move(-0.35, -4, -70.0);
    board2->RotateY(180);
    scene->AddSolidObject(board2);
    ConcreteBlock* block = new ConcreteBlock(Vector(15, 0.0, -50), Optics(Color(0.5, 0.5, 0.5)));
    block->RotateX(-0);
    block->RotateY(-15.0);
    block->RotateZ(-10);
    scene->AddSolidObject(block);
    ConcreteBlock* block2 = new ConcreteBlock(Vector(-15, 0.0, -50), Optics(Color(0.5, 0.5, 0.5)));
    block2->RotateX(-0);
    block2->RotateY(-15.0);
    block2->RotateZ(-100);
    scene->AddSolidObject(block2);
    for (size_t i = 0; i < 3; ++i)
    {
        Sphere* sphere = new Sphere(Vector(-1.8 + (1.8 * i), 1.5, -7), 0.5);
        sphere->SetOpacity(0.17 + (i * 0.3));
        sphere->SetMatteGlossBalance(0.95, Color(0.4, 0.5, 0.7), Color(0.8, 1.0, 0.7));
        scene->AddSolidObject(sphere);
    }
    Color co[10] = { Color(1, 1, 1),  Color(0, 1, 0), Color(1, 0, 0), Color(.4, .1, .1), Color(.5, .1, 1), Color(.1, .8, 1), Color(1, 1, .1), Color(.5, .5, .5), Color(1, .9, .9), Color(.5, .1, .5), };
    for (size_t i = 0; i < 10; ++i)
    {
        Sphere* c = new Sphere(Vector(0, 0, 0), .4);
        c->Move(Vector((-3.5 + (.9 * i)), -2, -7 - (7 * (int(i) % 2))));
        c->SetOpacity(1);
        c->SetMatteGlossBalance(0+(.1*i), co[i], co[i]);
        scene->AddSolidObject(c);
    }
    scene->AddLightSource(LightSource(Vector(+20.0, +20.0, +80.0), Color(0.5, 0.1, 0.1, 0.15)));
    scene->AddLightSource(LightSource(Vector(+100.0, +120.0, -70.0), Color(0.2, 0.5, 0.4, 1.00)));
    scene->AddLightSource(LightSource(Vector(+3.0, +13.0, +80.0), Color(0.6, 0.5, 0.3, 1.20)));
}
int main(int argc, char** argv)
{
    using namespace Imager;
    int procCount, rank; 
    std::chrono::system_clock::time_point start, end;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Imager::Scene scene(Color(0.37, 0.45, 0.37, 7.0e-6));
    Imager::ImageBuffer imagebuffer(WIDTH * AA_FACTOR, HEIGHT * AA_FACTOR, Color(0.37, 0.45, 0.37, 7.0e-6));
    constructLocalScene(&scene);
    MPI_Datatype PixelDataMpi, ColorMpi;
    int lengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
    MPI_Aint displacements[3];
    displacements[0] = offsetof(Imager::Color, red);
    displacements[1] = offsetof(Imager::Color, green);
    displacements[2] = offsetof(Imager::Color, blue);
    MPI_Type_create_struct(3, lengths, displacements, types, &ColorMpi);
    MPI_Type_commit(&ColorMpi);
    displacements[0] = offsetof(Imager::PixelData, color);
    displacements[1] = offsetof(Imager::PixelData, isAmbiguous);
    MPI_Datatype typess[2] = { ColorMpi, MPI_C_BOOL};
    MPI_Type_create_struct(2, lengths, displacements, typess, &PixelDataMpi);
    MPI_Type_commit(&PixelDataMpi);
    
    // Generate a PNG file that displays the scene...
    const char* filename = "RayTracedRender.png";
    const char* mt_filename = "Multi-Threaded RayTracedRender.png";
    const char* mpi_filename = "Distributed RayTracedRender.png";
    Imager::PixelData* GlobalBuffer;
    int imageArraySize = WIDTH*HEIGHT*AA_FACTOR*AA_FACTOR;
    if(rank==0) GlobalBuffer = new Imager::PixelData[imageArraySize];
    int lBufferSize = (imageArraySize)/procCount;
    Imager::PixelData* LocalBuffer = new Imager::PixelData[lBufferSize];
    printf("local buffer size: %d bytes\n", sizeof(Imager::PixelData)* lBufferSize);
    if(rank==0) start = std::chrono::high_resolution_clock::now();
    MPI_Scatter(GlobalBuffer, lBufferSize, PixelDataMpi, LocalBuffer, lBufferSize, PixelDataMpi, 0, MPI_COMM_WORLD);
    MPI_TraceRay(rank, procCount, LocalBuffer, &scene, WIDTH, HEIGHT, ZOOM, AA_FACTOR);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0) printf("Ray Trace calculation done");
    MPI_Gather(LocalBuffer, lBufferSize, PixelDataMpi, GlobalBuffer, lBufferSize, PixelDataMpi, 0, MPI_COMM_WORLD);
    if(rank==0)
    {
        imagebuffer.array = GlobalBuffer;
        SaveImageFromBuffer(mpi_filename, WIDTH, HEIGHT, ZOOM, AA_FACTOR, &imagebuffer);    
        end = std::chrono::high_resolution_clock::now();
        auto MPITime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();    
        printf("W: %d H: %d Anti-Aliase Factor: %d RecursionDepth: %d Tame take: %d milliseconds", 1920, 1080, AA_FACTOR, MT_MAX_OPTICAL_RECURSION_DEPTH, int(MPITime));

    }

    MPI_Finalize();
    return 1;
}

