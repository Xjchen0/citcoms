//==============================================================================
// citcoms_combine.cpp : Modern, high-performance C++17 postprocessing & combiner
// Replaces and accelerates legacy Python scripts:
//   - pasteCitcomData.py
//   - combine.py
//   - combinesurf.py
//   - autocombine.py / batchcombine.py
//
// Fully compatible with existing CitcomS data formats and parameters.
//==============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <algorithm>
#include <iomanip>

namespace fs = std::filesystem;

namespace CitcomsPost {

// Trim whitespace from string_view
inline std::string_view trim(std::string_view sv) {
    while (!sv.empty() && (sv.front() == ' ' || sv.front() == '\t' || sv.front() == '\r' || sv.front() == '\n')) {
        sv.remove_prefix(1);
    }
    while (!sv.empty() && (sv.back() == ' ' || sv.back() == '\t' || sv.back() == '\r' || sv.back() == '\n')) {
        sv.remove_suffix(1);
    }
    return sv;
}

// Split string_view by delimiter
inline std::vector<std::string> split(std::string_view sv, char delim) {
    std::vector<std::string> tokens;
    size_t start = 0;
    while (start < sv.size()) {
        size_t end = sv.find(delim, start);
        if (end == std::string_view::npos) {
            tokens.emplace_back(sv.substr(start));
            break;
        }
        tokens.emplace_back(sv.substr(start, end - start));
        start = end + 1;
    }
    return tokens;
}

// Get number of header lines to skip for a field infix
inline int get_header_lines(std::string_view infix) {
    if (infix == "coord" || infix == "botm" || infix == "surf" || infix == "visc") return 1;
    if (infix == "comp_nd" || infix == "pressure" || infix == "stress" || infix == "velo") return 2;
    return 0;
}

// Read lines from a file, skipping 'skip_headers' lines
std::vector<std::string> read_lines(const fs::path& filepath, int skip_headers) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath.string());
    }
    std::string line;
    for (int i = 0; i < skip_headers && std::getline(file, line); ++i) {
        // Skip header lines
    }
    std::vector<std::string> lines;
    lines.reserve(1024);
    while (std::getline(file, line)) {
        std::string_view sv = trim(line);
        if (!sv.empty()) {
            lines.emplace_back(sv);
        }
    }
    return lines;
}

// -----------------------------------------------------------------------------
// 1. Paste fields together for a single processor rank
// -----------------------------------------------------------------------------
bool paste_rank(const std::string& datadir, const std::string& datafile,
                const std::string& fields, int rank, int step,
                const std::string& savedir) {
    auto field_list = split(fields, ',');
    if (field_list.empty()) return false;

    std::vector<std::vector<std::string>> all_lines;
    all_lines.reserve(field_list.size());

    for (const auto& infix : field_list) {
        fs::path p;
        if (infix == "coord") {
            p = fs::path(datadir) / (datafile + ".coord." + std::to_string(rank));
        } else {
            p = fs::path(datadir) / (datafile + "." + infix + "." + std::to_string(rank) + "." + std::to_string(step));
        }
        int skip = get_header_lines(infix);
        try {
            all_lines.push_back(read_lines(p, skip));
        } catch (const std::exception& e) {
            std::cerr << "Warning: " << e.what() << std::endl;
            return false;
        }
    }

    size_t count = all_lines[0].size();
    for (size_t i = 1; i < all_lines.size(); ++i) {
        if (all_lines[i].size() != count) {
            std::cerr << "Error: Field line count mismatch for rank " << rank << std::endl;
            return false;
        }
    }

    fs::path outpath = fs::path(savedir) / (datafile + "." + fields + "." + std::to_string(rank) + "." + std::to_string(step) + ".pasted");
    std::ofstream out(outpath);
    if (!out.is_open()) {
        std::cerr << "Error: Cannot write to " << outpath << std::endl;
        return false;
    }

    // Set a large 64KB stream buffer for speed
    std::vector<char> buf(65536);
    out.rdbuf()->pubsetbuf(buf.data(), buf.size());

    for (size_t line_idx = 0; line_idx < count; ++line_idx) {
        for (size_t f = 0; f < all_lines.size(); ++f) {
            if (f > 0) out << ' ';
            out << all_lines[f][line_idx];
        }
        out << '\n';
    }
    return true;
}

// -----------------------------------------------------------------------------
// 2. Volume Combiner (replacing combine.py)
// -----------------------------------------------------------------------------
class VolumeCombiner {
public:
    int nodex, nodey, nodez;
    int nprocx, nprocy, nprocz;
    std::vector<std::string> saved;

    VolumeCombiner(int nx, int ny, int nz, int px, int py, int pz)
        : nodex(nx), nodey(ny), nodez(nz), nprocx(px), nprocy(py), nprocz(pz),
          saved(static_cast<size_t>(nx) * ny * nz) {}

    void join(const std::vector<std::string>& data, int me) {
        int mylocz = me % nprocz;
        int mylocx = ((me - mylocz) / nprocz) % nprocx;
        int mylocy = (((me - mylocz) / nprocz - mylocx) / nprocx) % nprocy;

        int mynodex = 1 + (nodex - 1) / nprocx;
        int mynodey = 1 + (nodey - 1) / nprocy;
        int mynodez = 1 + (nodez - 1) / nprocz;

        if (data.size() != static_cast<size_t>(mynodex * mynodey * mynodez)) {
            throw std::runtime_error("Incorrect data size for processor " + std::to_string(me));
        }

        int mynxs = (mynodex - 1) * mylocx;
        int mynys = (mynodey - 1) * mylocy;
        int mynzs = (mynodez - 1) * mylocz;

        size_t n = 0;
        for (int i = mynys; i < mynys + mynodey; ++i) {
            for (int j = mynxs; j < mynxs + mynodex; ++j) {
                for (int k = mynzs; k < mynzs + mynodez; ++k) {
                    size_t m = static_cast<size_t>(k) + j * nodez + i * nodex * nodez;
                    saved[m] = data[n++];
                }
            }
        }
    }

    void join_EL(const std::vector<std::string>& data, int me) {
        int mylocz = me % nprocz;
        int mylocx = ((me - mylocz) / nprocz) % nprocx;
        int mylocy = (((me - mylocz) / nprocz - mylocx) / nprocx) % nprocy;

        int mynodex = 1 + (nodex - nprocx) / nprocx;
        int mynodey = 1 + (nodey - nprocy) / nprocy;
        int mynodez = 1 + (nodez - nprocz) / nprocz;

        if (data.size() != static_cast<size_t>(mynodex * mynodey * mynodez)) {
            throw std::runtime_error("Incorrect EL data size for processor " + std::to_string(me));
        }

        int mynxs = mynodex * mylocx;
        int mynys = mynodey * mylocy;
        int mynzs = mynodez * mylocz;

        size_t n = 0;
        for (int i = mynys; i < mynys + mynodey; ++i) {
            for (int j = mynxs; j < mynxs + mynodex; ++j) {
                for (int k = mynzs; k < mynzs + mynodez; ++k) {
                    size_t m = static_cast<size_t>(k) + j * nodez + i * nodex * nodez;
                    saved[m] = data[n++];
                }
            }
        }
    }

    void write(const std::string& filename) const {
        std::ofstream out(filename);
        if (!out.is_open()) {
            throw std::runtime_error("Cannot write to output file: " + filename);
        }
        std::vector<char> buf(65536);
        out.rdbuf()->pubsetbuf(buf.data(), buf.size());

        out << nodex << " x " << nodey << " x " << nodez << '\n';
        for (const auto& line : saved) {
            out << line << '\n';
        }
    }
};

bool combine_volume(const std::string& prefix, int step,
                    int nodex, int nodey, int nodez,
                    int ncap, int nprocx, int nprocy, int nprocz,
                    const std::string& opts = "coord,velo,visc",
                    const std::string& datadir = ".") {
    int nproc_per_cap = nprocx * nprocy * nprocz;

    for (int i = 0; i < ncap; ++i) {
        VolumeCombiner cb(nodex, nodey, nodez, nprocx, nprocy, nprocz);
        for (int n = i * nproc_per_cap; n < (i + 1) * nproc_per_cap; ++n) {
            fs::path pasted = fs::path(datadir) / (prefix + "." + opts + "." + std::to_string(n) + "." + std::to_string(step) + ".pasted");
            // If pasted file does not exist, paste on the fly
            if (!fs::exists(pasted)) {
                paste_rank(datadir, prefix, opts, n, step, datadir);
            }
            std::cout << "reading " << pasted.string() << std::endl;
            auto data = read_lines(pasted, 0);
            if (opts == "heating") {
                cb.join_EL(data, n);
            } else {
                cb.join(data, n);
            }
        }

        std::ostringstream outname;
        if (opts == "coord,velo,visc") {
            outname << prefix << ".cap" << std::setfill('0') << std::setw(2) << i << "." << step;
        } else if (opts == "heating") {
            outname << prefix << ".ele" << std::setfill('0') << std::setw(2) << i << "." << step;
        } else {
            outname << prefix << ".opt" << std::setfill('0') << std::setw(2) << i << "." << step;
        }

        fs::path outpath = fs::path(datadir) / outname.str();
        std::cout << "writing " << outpath.string() << std::endl;
        cb.write(outpath.string());
    }
    return true;
}

// -----------------------------------------------------------------------------
// 3. Surface Combiner (replacing combinesurf.py)
// -----------------------------------------------------------------------------
bool combine_surface(const std::string& prefix, int step,
                     int nox, int noy, int noz,
                     int nprocxy, int nprocx, int nprocy, int nprocz,
                     const std::string& datadir = ".") {
    int nproc_per_cap = nprocx * nprocy * nprocz;
    int mynox = 1 + (nox - 1) / nprocx;
    int mynoy = 1 + (noy - 1) / nprocy;
    int mynoz = 1 + (noz - 1) / nprocz;
    int snodes = mynox * mynoy;

    for (int i = 0; i < nprocxy; ++i) {
        std::vector<std::string> saved(static_cast<size_t>(nox) * noy);

        for (int n = i * nproc_per_cap + (nprocz - 1); n < (i + 1) * nproc_per_cap; n += nprocz) {
            fs::path crdfile = fs::path(datadir) / (prefix + ".coord." + std::to_string(n));
            fs::path surffile = fs::path(datadir) / (prefix + ".surf." + std::to_string(n) + "." + std::to_string(step));

            std::cout << "reading " << surffile.string() << std::endl;
            auto crd = read_lines(crdfile, 1);
            auto surf = read_lines(surffile, 1);

            int mylocz = n % nprocz;
            int mylocx = ((n - mylocz) / nprocz) % nprocx;
            int mylocy = (((n - mylocz) / nprocz - mylocx) / nprocx) % nprocy;

            int mynxs = (mynox - 1) * mylocx;
            int mynys = (mynoy - 1) * mylocy;

            size_t local_idx = 0;
            for (int iy = mynys; iy < mynys + mynoy; ++iy) {
                for (int ix = mynxs; ix < mynxs + mynox; ++ix) {
                    size_t m = ix + static_cast<size_t>(iy) * nox;
                    // Coordinate line at surface
                    size_t crd_line = (local_idx + 1) * mynoz - 1;
                    if (crd_line < crd.size() && local_idx < surf.size()) {
                        // Extract x, y
                        auto tokens = split(crd[crd_line], ' ');
                        std::string x = tokens.size() > 0 ? tokens[0] : "0.0";
                        std::string y = tokens.size() > 1 ? tokens[1] : "0.0";
                        saved[m] = x + " " + y + " " + surf[local_idx];
                    }
                    local_idx++;
                }
            }
        }

        std::string outname = prefix + ".surf" + std::to_string(i) + "." + std::to_string(step);
        fs::path outpath = fs::path(datadir) / outname;
        std::cout << "writing " << outpath.string() << std::endl;

        std::ofstream out(outpath);
        out << nox << " x " << noy << '\n';
        for (const auto& line : saved) {
            out << line << '\n';
        }
    }
    return true;
}

} // namespace CitcomsPost

void print_usage(const char* prog) {
    std::cout << "CitcomS Modern C++17 Combiner & Postprocessor\n"
              << "Usage:\n"
              << "  1. Volume Combine (replaces combine.py):\n"
              << "     " << prog << " volume <prefix> <step> <nodex> <nodey> <nodez> <ncap> <nprocx> <nprocy> <nprocz> [opts=coord,velo,visc] [datadir=.]\n\n"
              << "  2. Surface Combine (replaces combinesurf.py):\n"
              << "     " << prog << " surf <prefix> <step> <nodex> <nodey> <nodez> <nprocxy> <nprocx> <nprocy> <nprocz> [datadir=.]\n\n"
              << "  3. Paste Fields (replaces pasteCitcomData.py):\n"
              << "     " << prog << " paste <datadir> <datafile> <fields> <rank> <step> <savedir>\n"
              << std::endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    std::string mode = argv[1];

    try {
        if (mode == "volume") {
            if (argc < 11) {
                print_usage(argv[0]);
                return 1;
            }
            std::string prefix = argv[2];
            int step = std::stoi(argv[3]);
            int nodex = std::stoi(argv[4]);
            int nodey = std::stoi(argv[5]);
            int nodez = std::stoi(argv[6]);
            int ncap = std::stoi(argv[7]);
            int nprocx = std::stoi(argv[8]);
            int nprocy = std::stoi(argv[9]);
            int nprocz = std::stoi(argv[10]);
            std::string opts = (argc >= 12) ? argv[11] : "coord,velo,visc";
            std::string datadir = (argc >= 13) ? argv[12] : ".";

            CitcomsPost::combine_volume(prefix, step, nodex, nodey, nodez, ncap, nprocx, nprocy, nprocz, opts, datadir);
            return 0;
        } else if (mode == "surf") {
            if (argc < 11) {
                print_usage(argv[0]);
                return 1;
            }
            std::string prefix = argv[2];
            int step = std::stoi(argv[3]);
            int nox = std::stoi(argv[4]);
            int noy = std::stoi(argv[5]);
            int noz = std::stoi(argv[6]);
            int nprocxy = std::stoi(argv[7]);
            int nprocx = std::stoi(argv[8]);
            int nprocy = std::stoi(argv[9]);
            int nprocz = std::stoi(argv[10]);
            std::string datadir = (argc >= 12) ? argv[11] : ".";

            CitcomsPost::combine_surface(prefix, step, nox, noy, noz, nprocxy, nprocx, nprocy, nprocz, datadir);
            return 0;
        } else if (mode == "paste") {
            if (argc < 8) {
                print_usage(argv[0]);
                return 1;
            }
            std::string datadir = argv[2];
            std::string datafile = argv[3];
            std::string fields = argv[4];
            int rank = std::stoi(argv[5]);
            int step = std::stoi(argv[6]);
            std::string savedir = argv[7];

            bool ok = CitcomsPost::paste_rank(datadir, datafile, fields, rank, step, savedir);
            return ok ? 0 : 1;
        } else {
            // Also support direct combine.py argument syntax:
            // <prefix> <step> <nodex> <nodey> <nodez> <ncap> <nprocx> <nprocy> <nprocz>
            if (argc == 10) {
                std::string prefix = argv[1];
                int step = std::stoi(argv[2]);
                int nodex = std::stoi(argv[3]);
                int nodey = std::stoi(argv[4]);
                int nodez = std::stoi(argv[5]);
                int ncap = std::stoi(argv[6]);
                int nprocx = std::stoi(argv[7]);
                int nprocy = std::stoi(argv[8]);
                int nprocz = std::stoi(argv[9]);
                CitcomsPost::combine_volume(prefix, step, nodex, nodey, nodez, ncap, nprocx, nprocy, nprocz);
                return 0;
            }
            std::cerr << "Unknown mode: " << mode << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
