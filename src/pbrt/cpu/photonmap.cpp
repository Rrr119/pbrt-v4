// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include <pbrt/cpu/photonmap.h>

#include <pbrt/bsdf.h>
#include <pbrt/cameras.h>
#include <pbrt/film.h>
#include <pbrt/interaction.h>
#include <pbrt/lights.h>
#include <pbrt/lightsamplers.h>
#include <pbrt/options.h>
#include <pbrt/paramdict.h>
#include <pbrt/samplers.h>
#include <pbrt/util/check.h>
#include <pbrt/util/color.h>
#include <pbrt/util/colorspace.h>
#include <pbrt/util/display.h>
#include <pbrt/util/error.h>
#include <pbrt/util/hash.h>
#include <pbrt/util/math.h>
#include <pbrt/util/memory.h>
#include <pbrt/util/parallel.h>
#include <pbrt/util/print.h>
#include <pbrt/util/progressreporter.h>
#include <pbrt/util/pstd.h>
#include <pbrt/util/rng.h>
#include <pbrt/util/spectrum.h>
#include <pbrt/util/string.h>

#include <algorithm>
#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

namespace pbrt {

STAT_COUNTER("Photon map/Global estimate calls", nGlobalPhotonEstimates);
STAT_COUNTER("Photon map/Global photons found", nGlobalPhotonsFound);
STAT_COUNTER("Photon map/Global nonzero estimates", nGlobalNonzeroEstimates);
STAT_COUNTER("Photon map/Caustic estimate calls", nCausticPhotonEstimates);
STAT_COUNTER("Photon map/Caustic photons found", nCausticPhotonsFound);
STAT_COUNTER("Photon map/Caustic nonzero estimates", nCausticNonzeroEstimates);

// Photon Definition
struct Photon {
    Point3f p;
    Vector3f wi;
    RGB power;
};

struct PhotonCellHash {
    size_t operator()(Point3i p) const {
        return size_t(Hash(p.x, p.y, p.z));
    }
};

struct PhotonCandidate {
    Float distanceSquared;
    const Photon *photon;
};

// PhotonMapData Definition
class PhotonMapData {
  public:
    explicit PhotonMapData(bool caustic) : caustic(caustic) {}

    void Build(std::vector<Photon> photons, Float radius);

    RGB Estimate(const Point3f &p, const Normal3f &n, const Vector3f &wo,
                 const BSDF &bsdf, const SampledWavelengths &lambda,
                 const RGBColorSpace *colorSpace, Float radius, int maxPhotons) const;

  private:
    // PhotonMapData Private Members
    std::vector<Photon> photons;
    std::unordered_map<Point3i, std::vector<uint32_t>, PhotonCellHash> grid;
    Float cellSize = 1;
    bool caustic;
};

void PhotonMapData::Build(std::vector<Photon> p, Float radius) {
    photons = std::move(p);
    grid.clear();
    cellSize = std::max<Float>(radius, 1e-5f);

    grid.reserve(photons.size());
    Bounds3f bounds;
    RGB powerSum(0, 0, 0);
    for (uint32_t i = 0; i < photons.size(); ++i) {
        const Point3f &photonPosition = photons[i].p;
        bounds = Union(bounds, photonPosition);
        powerSum += photons[i].power;
        Point3i cell(pstd::floor(photonPosition.x / cellSize),
                     pstd::floor(photonPosition.y / cellSize),
                     pstd::floor(photonPosition.z / cellSize));
        grid[cell].push_back(i);
    }
    LOG_VERBOSE("%s photon map: %d photons, bounds %s, average power %s",
                caustic ? "Caustic" : "Global", int(photons.size()), bounds,
                photons.empty() ? RGB(0, 0, 0) : powerSum / Float(photons.size()));
}

RGB PhotonMapData::Estimate(const Point3f &p, const Normal3f &n, const Vector3f &wo,
                            const BSDF &bsdf, const SampledWavelengths &lambda,
                            const RGBColorSpace *colorSpace, Float radius,
                            int maxPhotons) const {
    if (caustic)
        ++nCausticPhotonEstimates;
    else
        ++nGlobalPhotonEstimates;

    if (photons.empty() || radius <= 0)
        return RGB(0, 0, 0);

    Point3i centerCell(pstd::floor(p.x / cellSize), pstd::floor(p.y / cellSize),
                       pstd::floor(p.z / cellSize));
    Float radiusSquared = Sqr(radius);

    std::vector<PhotonCandidate> candidates;
    candidates.reserve(64);
    for (int dz = -1; dz <= 1; ++dz) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                Point3i cell(centerCell.x + dx, centerCell.y + dy, centerCell.z + dz);
                auto iter = grid.find(cell);
                if (iter == grid.end())
                    continue;
                for (uint32_t photonIndex : iter->second) {
                    const Photon &photon = photons[photonIndex];
                    Float distanceSquared = DistanceSquared(p, photon.p);
                    if (distanceSquared <= radiusSquared)
                        candidates.push_back({distanceSquared, &photon});
                }
            }
        }
    }
    if (candidates.empty())
        return RGB(0, 0, 0);
    if (caustic)
        nCausticPhotonsFound += candidates.size();
    else
        nGlobalPhotonsFound += candidates.size();

    // Use the nearest photons and the distance to the last one as the
    // effective kernel radius, as in the standard photon-map estimate.
    int nPhotons = int(candidates.size());
    if (maxPhotons > 0 && nPhotons > maxPhotons) {
        std::nth_element(candidates.begin(), candidates.begin() + maxPhotons,
                         candidates.end(), [](const PhotonCandidate &a,
                                              const PhotonCandidate &b) {
                             return a.distanceSquared < b.distanceSquared;
                         });
        nPhotons = maxPhotons;
        radiusSquared = candidates[maxPhotons - 1].distanceSquared;
    }
    if (radiusSquared <= 0)
        return RGB(0, 0, 0);

    RGB sum(0, 0, 0);
    for (int i = 0; i < nPhotons; ++i) {
        const Photon &photon = *candidates[i].photon;
        SampledSpectrum f = bsdf.f(wo, photon.wi) * AbsDot(photon.wi, n);
        if (!f)
            continue;
        RGB fRGB = f.ToRGB(lambda, *colorSpace);
        sum += photon.power * fRGB;
    }
    RGB estimate = sum / (Pi * radiusSquared);
    if (estimate != RGB(0, 0, 0)) {
        if (caustic)
            ++nCausticNonzeroEstimates;
        else
            ++nGlobalNonzeroEstimates;
    }
    return estimate;
}

namespace {

// Photon Shooting Helper
std::vector<Photon> ShootPhotons(const Integrator &integrator, Camera camera,
                                 Sampler samplerPrototype, int maxPhotonDepth,
                                 int nPhotons, bool causticMap, const Film &film,
                                 const RGBColorSpace *colorSpace,
                                 PowerLightSampler &lightSampler, int seed,
                                 ProgressReporter &progress) {
    if (nPhotons <= 0)
        return {};

    ThreadLocal<std::vector<Photon>> photonLists(
        []() { return std::vector<Photon>(); });
    ThreadLocal<ScratchBuffer> scratchBuffers(
        []() { return ScratchBuffer(1024 * 1024); });
    ThreadLocal<Sampler> samplers(
        [&]() { return samplerPrototype.Clone(Allocator()); });

    ParallelFor(0, nPhotons, [&](int64_t start, int64_t end) {
        std::vector<Photon> &localPhotons = photonLists.Get();
        ScratchBuffer &scratchBuffer = scratchBuffers.Get();
        Sampler sampler = samplers.Get();
        Float invPhotonPaths = 1.f / Float(nPhotons);

        for (int64_t photonIndex = start; photonIndex < end; ++photonIndex) {
            sampler.StartPixelSample(Point2i(0, 0),
                                     int(photonIndex & 0x7fffffff));
            RNG rng(Hash(uint64_t(photonIndex), uint64_t(seed),
                         uint64_t(causticMap)));

            // Choose a wavelength and a light source for this photon.
            Float uLambda = Options->disableWavelengthJitter
                                ? Float(0.5)
                                : rng.Uniform<Float>();
            SampledWavelengths lambda = film.SampleWavelengths(uLambda);

            pstd::optional<SampledLight> sampledLight =
                lightSampler.Sample(rng.Uniform<Float>());
            if (!sampledLight)
                continue;
            Light light = sampledLight->light;
            Float p_l = sampledLight->p;

            Float time = camera.SampleTime(rng.Uniform<Float>());
            pstd::optional<LightLeSample> les =
                light.SampleLe({rng.Uniform<Float>(), rng.Uniform<Float>()},
                               {rng.Uniform<Float>(), rng.Uniform<Float>()}, lambda,
                               time);
            if (!les || les->pdfPos == 0 || les->pdfDir == 0 || !les->L)
                continue;

            RayDifferential photonRay(les->ray);
            SampledSpectrum beta = les->L * les->AbsCosTheta(photonRay.d) /
                                   (p_l * les->pdfPos * les->pdfDir);
            if (!beta)
                continue;

            bool hasSpecular = false;
            for (int depth = 0; depth < maxPhotonDepth; ++depth) {
                pstd::optional<ShapeIntersection> si = integrator.Intersect(photonRay);
                if (!si)
                    break;
                SurfaceInteraction &isect = si->intr;

                BSDF bsdf = isect.GetBSDF(photonRay, lambda, camera, scratchBuffer,
                                          sampler);
                if (!bsdf) {
                    isect.SkipIntersection(&photonRay, si->tHit);
                    --depth;
                    continue;
                }

                BxDFFlags flags = bsdf.Flags();
                if (IsDiffuse(flags)) {
                    // Do not store direct photons; direct illumination is
                    // sampled explicitly along camera paths.
                    if (depth > 0) {
                        if (causticMap) {
                            if (hasSpecular) {
                                localPhotons.push_back(
                                    {isect.p(), -photonRay.d,
                                     (beta * invPhotonPaths)
                                         .ToRGB(lambda, *colorSpace)});
                                break;
                            }
                        } else if (!hasSpecular)
                            localPhotons.push_back(
                                {isect.p(), -photonRay.d,
                                 (beta * invPhotonPaths).ToRGB(lambda, *colorSpace)});
                    }

                    // Caustic photons are only stored at the first diffuse
                    // surface after a specular event.  A caustic path that
                    // later reaches another diffuse surface does not belong
                    // to the global photon map.
                    if (hasSpecular)
                        break;
                }

                // Sample a new photon direction.
                Vector3f wo = -photonRay.d;
                pstd::optional<BSDFSample> bs = bsdf.Sample_f(
                    wo, rng.Uniform<Float>(),
                    {rng.Uniform<Float>(), rng.Uniform<Float>()},
                    TransportMode::Importance);
                if (!bs)
                    break;
                hasSpecular |= bs->IsSpecular();

                SampledSpectrum bnew =
                    beta * bs->f * AbsDot(bs->wi, isect.shading.n) / bs->pdf;
                if (!bnew)
                    break;

                // Russian roulette on photon throughput.
                Float q = std::max<Float>(
                    0, 1 - bnew.MaxComponentValue() / beta.MaxComponentValue());
                if (rng.Uniform<Float>() < q)
                    break;
                beta = bnew / (1 - q);
                photonRay = RayDifferential(isect.SpawnRay(bs->wi));
            }
        }
        progress.Update(end - start);
        scratchBuffer.Reset();
    });

    std::vector<Photon> photons;
    photonLists.ForAll([&](std::vector<Photon> &localPhotons) {
        photons.insert(photons.end(), localPhotons.begin(), localPhotons.end());
    });
    return photons;
}

}  // namespace

// PhotonMapIntegrator Method Definitions
std::unique_ptr<PhotonMapIntegrator> PhotonMapIntegrator::Create(
    const ParameterDictionary &parameters, const RGBColorSpace *colorSpace, Camera camera,
    Sampler sampler, Primitive aggregate, std::vector<Light> lights, const FileLoc *loc) {
    int maxDepth = parameters.GetOneInt("maxdepth", 5);
    int maxPhotonDepth = parameters.GetOneInt("maxphotondepth", 10);

    int globalPhotonCount = parameters.GetOneInt(
        "photons", parameters.GetOneInt("globalphotons", 100000));
    int causticPhotonCount = parameters.GetOneInt("causticphotons", 10000);

    Float globalRadius =
        parameters.GetOneFloat("radius", parameters.GetOneFloat("maxdist", 0.1f));
    Float causticRadius = parameters.GetOneFloat(
        "causticradius", parameters.GetOneFloat("causticmaxdist", globalRadius));
    int maxPhotons = parameters.GetOneInt("maxphotons", 50);
    bool caustics = parameters.GetOneBool("caustics", true);
    int seed = parameters.GetOneInt("seed", Options->seed);

    if (maxDepth < 0)
        ErrorExit(loc, "\"maxdepth\" must be non-negative for the photon-map integrator.");
    if (maxPhotonDepth < 0)
        ErrorExit(loc, "\"maxphotondepth\" must be non-negative for the photon-map "
                       "integrator.");
    if (globalPhotonCount < 0 || causticPhotonCount < 0)
        ErrorExit(loc, "Photon counts must be non-negative.");
    if (globalRadius <= 0 || causticRadius <= 0)
        ErrorExit(loc, "Photon-map search radii must be positive.");
    if (maxPhotons < 0)
        ErrorExit(loc, "\"maxphotons\" must be non-negative.");

    if (!camera.GetFilm().Is<RGBFilm>() && !camera.GetFilm().Is<GBufferFilm>())
        ErrorExit(loc, "The photon-map integrator requires an RGB or gbuffer film.");

    return std::make_unique<PhotonMapIntegrator>(
        camera, sampler, aggregate, lights, maxDepth, maxPhotonDepth,
        globalPhotonCount, causticPhotonCount, globalRadius, causticRadius, maxPhotons,
        caustics, seed, colorSpace);
}

std::string PhotonMapIntegrator::ToString() const {
    return StringPrintf(
        "[ PhotonMapIntegrator camera: %s maxDepth: %d maxPhotonDepth: %d "
        "globalPhotons: %d causticPhotons: %d radius: %f causticRadius: %f "
        "maxPhotons: %d caustics: %s seed: %d colorSpace: %s ]",
        camera, maxDepth, maxPhotonDepth, globalPhotonCount, causticPhotonCount,
        globalRadius, causticRadius, maxPhotons, caustics, seed, *colorSpace);
}

SampledSpectrum PhotonMapIntegrator::SampleLd(const SurfaceInteraction &intr,
                                              const BSDF &b, SampledWavelengths &lambda,
                                              Sampler sampler,
                                              LightSampler lightSampler) const {
    const BSDF *bsdf = &b;
    LightSampleContext ctx(intr);
    BxDFFlags flags = bsdf->Flags();
    if (IsReflective(flags) && !IsTransmissive(flags))
        ctx.pi = intr.OffsetRayOrigin(intr.wo);
    else if (IsTransmissive(flags) && !IsReflective(flags))
        ctx.pi = intr.OffsetRayOrigin(-intr.wo);

    pstd::optional<SampledLight> sampledLight =
        lightSampler.Sample(ctx, sampler.Get1D());
    if (!sampledLight)
        return {};

    Light light = sampledLight->light;
    DCHECK(light && sampledLight->p > 0);
    pstd::optional<LightLiSample> ls =
        light.SampleLi(ctx, sampler.Get2D(), lambda, true);
    if (!ls || !ls->L || ls->pdf == 0)
        return {};

    Vector3f wo = intr.wo, wi = ls->wi;
    SampledSpectrum f = bsdf->f(wo, wi) * AbsDot(wi, intr.shading.n);
    if (!f || !Unoccluded(intr, ls->pLight))
        return {};

    Float p_l = sampledLight->p * ls->pdf;
    if (IsDeltaLight(light.Type()))
        return ls->L * f / p_l;

    Float p_b = bsdf->PDF(wo, wi);
    Float w_l = PowerHeuristic(1, p_l, 1, p_b);
    return w_l * ls->L * f / p_l;
}

SampledSpectrum PhotonMapIntegrator::Li(RayDifferential ray,
                                        SampledWavelengths &lambda, Sampler sampler,
                                        ScratchBuffer &scratchBuffer,
                                        const PhotonMapData &globalMap,
                                        const PhotonMapData &causticMap,
                                        LightSampler lightSampler,
                                        VisibleSurface *visibleSurface) const {
    SampledSpectrum L(0.f), beta(1.f);
    int depth = 0;
    bool specularBounce = true;
    Float p_b = 0;
    LightSampleContext prevIntrCtx;

    while (true) {
        pstd::optional<ShapeIntersection> si = Intersect(ray);
        if (!si) {
            for (const Light &light : infiniteLights) {
                SampledSpectrum Le = light.Le(ray, lambda);
                if (depth == 0 || specularBounce)
                    L += beta * Le;
                else {
                    Float p_l = lightSampler.PMF(prevIntrCtx, light) *
                                light.PDF_Li(prevIntrCtx, ray.d, true);
                    Float w_b = PowerHeuristic(1, p_b, 1, p_l);
                    L += beta * w_b * Le;
                }
            }
            break;
        }

        SurfaceInteraction &isect = si->intr;
        SampledSpectrum Le = isect.Le(-ray.d, lambda);
        if (Le) {
            if (depth == 0 || specularBounce)
                L += beta * Le;
            else {
                Light areaLight(isect.areaLight);
                Float p_l = lightSampler.PMF(prevIntrCtx, areaLight) *
                            areaLight.PDF_Li(prevIntrCtx, ray.d, true);
                Float w_l = PowerHeuristic(1, p_b, 1, p_l);
                L += beta * w_l * Le;
            }
        }

        BSDF bsdf = isect.GetBSDF(ray, lambda, camera, scratchBuffer, sampler);
        if (!bsdf) {
            specularBounce = true;
            isect.SkipIntersection(&ray, si->tHit);
            continue;
        }

        Vector3f wo = -ray.d;
        if (depth == 0 && visibleSurface) {
            constexpr int nRhoSamples = 16;
            const Float ucRho[nRhoSamples] = {
                0.75741637, 0.37870818, 0.7083487,  0.18935409, 0.9149363,  0.35417435,
                0.5990858,  0.09467703, 0.8578725,  0.45746812, 0.686759,   0.17708716,
                0.9674518,  0.2995429,  0.5083201,  0.047338516};
            const Point2f uRho[nRhoSamples] = {
                Point2f(0.855985, 0.570367),  Point2f(0.381823, 0.851844),
                Point2f(0.285328, 0.764262),  Point2f(0.733380, 0.114073),
                Point2f(0.542663, 0.344465),  Point2f(0.127274, 0.414848),
                Point2f(0.964700, 0.947162),  Point2f(0.594089, 0.643463),
                Point2f(0.095109, 0.170369),  Point2f(0.825444, 0.263359),
                Point2f(0.429467, 0.454469),  Point2f(0.244460, 0.816459),
                Point2f(0.756135, 0.731258),  Point2f(0.516165, 0.152852),
                Point2f(0.180888, 0.214174),  Point2f(0.898579, 0.503897)};
            SampledSpectrum albedo = bsdf.rho(wo, ucRho, uRho);
            *visibleSurface = VisibleSurface(isect, albedo, lambda);
        }

        if (depth++ == maxDepth)
            break;

        // Direct illumination is sampled explicitly at every non-specular
        // camera-path vertex.
        if (IsNonSpecular(bsdf.Flags()))
            L += beta *
                 SampleLd(isect, bsdf, lambda, sampler, lightSampler);

        BxDFFlags flags = bsdf.Flags();
        if (IsDiffuse(flags)) {
            RGB photonEstimate(0, 0, 0);
            photonEstimate += globalMap.Estimate(
                isect.p(), isect.shading.n, wo, bsdf, lambda, colorSpace,
                globalRadius, maxPhotons);
            if (caustics)
                photonEstimate += causticMap.Estimate(
                    isect.p(), isect.shading.n, wo, bsdf, lambda, colorSpace,
                    causticRadius, maxPhotons);

            photonEstimate = ClampZero(photonEstimate);
            if ((photonEstimate.r != 0 || photonEstimate.g != 0 ||
                 photonEstimate.b != 0) &&
                IsFinite(photonEstimate.r) && IsFinite(photonEstimate.g) &&
                IsFinite(photonEstimate.b)) {
                SampledSpectrum spectrum =
                    RGBUnboundedSpectrum(*colorSpace, photonEstimate).Sample(lambda);
                L += beta * spectrum;
            }
            // The photon maps account for the remaining diffuse transport;
            // terminating here avoids double counting it with BSDF sampling.
            break;
        }

        // Glossy surfaces cannot use the diffuse photon estimate, but they
        // can still continue along the camera path.
        pstd::optional<BSDFSample> bs = bsdf.Sample_f(
            wo, sampler.Get1D(), sampler.Get2D(), TransportMode::Radiance);
        if (!bs)
            break;

        beta *= bs->f * AbsDot(bs->wi, isect.shading.n) / bs->pdf;
        if (!beta)
            break;
        p_b = bs->pdfIsProportional ? bsdf.PDF(wo, bs->wi) : bs->pdf;
        specularBounce = bs->IsSpecular();
        prevIntrCtx = LightSampleContext(isect);

        SampledSpectrum rrBeta = beta;
        if (rrBeta.MaxComponentValue() < 1) {
            Float q = std::max<Float>(0.05f, 1 - rrBeta.MaxComponentValue());
            if (sampler.Get1D() < q)
                break;
            beta /= 1 - q;
        }

        ray = isect.SpawnRay(ray, bsdf, bs->wi, bs->flags, bs->eta);
    }
    return L;
}

void PhotonMapIntegrator::Render() {
    Film film = camera.GetFilm();
    Bounds2i pixelBounds = film.PixelBounds();
    int spp = samplerPrototype.SamplesPerPixel();

    // Build the global and caustic photon maps before rendering camera paths.
    ProgressReporter photonProgress(globalPhotonCount + causticPhotonCount,
                                    "Photon mapping", Options->quiet);
    PowerLightSampler shootLightSampler(lights, Allocator());

    std::vector<Photon> globalPhotons = ShootPhotons(
        *this, camera, samplerPrototype, maxPhotonDepth, globalPhotonCount, false,
        film, colorSpace, shootLightSampler, seed, photonProgress);
    std::vector<Photon> causticPhotons;
    if (caustics && causticPhotonCount > 0)
        causticPhotons = ShootPhotons(
            *this, camera, samplerPrototype, maxPhotonDepth, causticPhotonCount, true,
            film, colorSpace, shootLightSampler, seed ^ 0x9e3779b9, photonProgress);
    else
        photonProgress.Update(causticPhotonCount);
    photonProgress.Done();

    LOG_VERBOSE("Photon map built with %d global and %d caustic photons",
                int(globalPhotons.size()), int(causticPhotons.size()));

    PhotonMapData globalMap(false), causticMap(true);
    globalMap.Build(std::move(globalPhotons), globalRadius);
    causticMap.Build(std::move(causticPhotons), causticRadius);

    BVHLightSampler lightSampler(lights, Allocator());
    ThreadLocal<ScratchBuffer> scratchBuffers([]() { return ScratchBuffer(); });
    ThreadLocal<Sampler> samplers(
        [&]() { return samplerPrototype.Clone(Allocator()); });

    ProgressReporter progress(int64_t(spp) * pixelBounds.Area(), "Rendering",
                              Options->quiet);
    ParallelFor2D(pixelBounds, [&](Bounds2i tileBounds) {
        ScratchBuffer &scratchBuffer = scratchBuffers.Get();
        Sampler sampler = samplers.Get();
        for (Point2i pPixel : tileBounds) {
            for (int sampleIndex = 0; sampleIndex < spp; ++sampleIndex) {
                sampler.StartPixelSample(pPixel, sampleIndex);

                Float lu = sampler.Get1D();
                if (Options->disableWavelengthJitter)
                    lu = 0.5f;
                SampledWavelengths lambda = film.SampleWavelengths(lu);

                CameraSample cameraSample =
                    GetCameraSample(sampler, pPixel, film.GetFilter());
                pstd::optional<CameraRayDifferential> cameraRay =
                    camera.GenerateRayDifferential(cameraSample, lambda);

                SampledSpectrum L(0.f);
                VisibleSurface visibleSurface;
                if (cameraRay) {
                    Float rayDiffScale =
                        std::max<Float>(.125f,
                                        1 / std::sqrt(Float(sampler.SamplesPerPixel())));
                    if (!Options->disablePixelJitter)
                        cameraRay->ray.ScaleDifferentials(rayDiffScale);

                    L = cameraRay->weight *
                        Li(cameraRay->ray, lambda, sampler, scratchBuffer, globalMap,
                           causticMap, &lightSampler,
                           film.UsesVisibleSurface() ? &visibleSurface : nullptr);
                }

                if (L.HasNaNs() || IsInf(L.y(lambda))) {
                    LOG_ERROR("Invalid radiance value returned for pixel (%d, %d), "
                              "sample %d. Setting to black.",
                              pPixel.x, pPixel.y, sampleIndex);
                    L = SampledSpectrum(0.f);
                }

                film.AddSample(pPixel, L, lambda, &visibleSurface,
                               cameraSample.filterWeight);
                scratchBuffer.Reset();
            }
            progress.Update(spp);
        }
    });
    progress.Done();

    ImageMetadata metadata;
    metadata.renderTimeSeconds = progress.ElapsedSeconds();
    metadata.samplesPerPixel = spp;
    metadata.pixelBounds = pixelBounds;
    metadata.fullResolution = film.FullResolution();
    metadata.colorSpace = colorSpace;
    camera.InitMetadata(&metadata);
    film.WriteImage(metadata, 1.f / spp);

    DisconnectFromDisplayServer();
}

}  // namespace pbrt
