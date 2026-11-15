// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_CPU_PHOTONMAP_H
#define PBRT_CPU_PHOTONMAP_H

#include <pbrt/cpu/integrators.h>

#include <pbrt/film.h>
#include <pbrt/util/progressreporter.h>

#include <memory>
#include <string>
#include <vector>

namespace pbrt {

class PhotonMapData;

// PhotonMapIntegrator Definition
class PhotonMapIntegrator : public Integrator {
  public:
    // PhotonMapIntegrator Public Methods
    PhotonMapIntegrator(Camera camera, Sampler sampler, Primitive aggregate,
                        std::vector<Light> lights, int maxDepth, int maxPhotonDepth,
                        int globalPhotonCount, int causticPhotonCount, Float globalRadius,
                        Float causticRadius, int maxPhotons, bool caustics, int seed,
                        const RGBColorSpace *colorSpace)
        : Integrator(aggregate, lights),
          camera(camera),
          samplerPrototype(sampler),
          maxDepth(maxDepth),
          maxPhotonDepth(maxPhotonDepth),
          globalPhotonCount(globalPhotonCount),
          causticPhotonCount(causticPhotonCount),
          globalRadius(globalRadius),
          causticRadius(causticRadius),
          maxPhotons(maxPhotons),
          caustics(caustics),
          seed(seed),
          colorSpace(colorSpace) {}

    static std::unique_ptr<PhotonMapIntegrator> Create(
        const ParameterDictionary &parameters, const RGBColorSpace *colorSpace,
        Camera camera, Sampler sampler, Primitive aggregate, std::vector<Light> lights,
        const FileLoc *loc);

    std::string ToString() const;

    void Render();

  private:
    // PhotonMapIntegrator Private Methods
    SampledSpectrum SampleLd(const SurfaceInteraction &intr, const BSDF &bsdf,
                             SampledWavelengths &lambda, Sampler sampler,
                             LightSampler lightSampler) const;

    SampledSpectrum Li(RayDifferential ray, SampledWavelengths &lambda, Sampler sampler,
                       ScratchBuffer &scratchBuffer, const PhotonMapData &globalMap,
                       const PhotonMapData &causticMap, LightSampler lightSampler,
                       VisibleSurface *visibleSurface) const;

    // PhotonMapIntegrator Private Members
    Camera camera;
    Sampler samplerPrototype;
    int maxDepth, maxPhotonDepth;
    int globalPhotonCount, causticPhotonCount;
    Float globalRadius, causticRadius;
    int maxPhotons;
    bool caustics;
    int seed;
    const RGBColorSpace *colorSpace;
};

}  // namespace pbrt

#endif  // PBRT_CPU_PHOTONMAP_H
