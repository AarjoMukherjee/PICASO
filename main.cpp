#include "Common.h"
#include "Calibration.h"
#include "CoincidenceProcessor.h"
#include "OutputManager.h"
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <map>
#include <sstream>
#include <chrono>
#include <utility>
#include <vector>
#include <string>
#include <ctime>

int extract_run_number(const std::string& fname)
{
    size_t modPos = fname.rfind("_mod");
    if (modPos == std::string::npos)
        return -1;

    size_t prevUnderscore = fname.rfind('_', modPos - 1);
    if (prevUnderscore == std::string::npos)
        return -1;

    return std::stoi(
        fname.substr(
            prevUnderscore + 1,
            modPos - prevUnderscore - 1
        )
    );
}


int main()
{
    auto t_start =
    std::chrono::high_resolution_clock::now();
    // int seed = time(NULL);
    // srand(seed);
    // printf("Random number seed: %d", seed);
    int Nfiles;

    double timewin;

    char file_list[MAX_FILES][150];

    // ========================================================
    // DEFINE CODE PARAMETERS
    // ========================================================

    printf("Parameter Values:");

    printf("\nTimeBinWidthPs: ");
    scanf("%f",&TimeBinWidthPs);

    printf("\nTimeOffsetBin: ");
    scanf("%f",&TimeOffsetBin);

    printf("\nEnergyScale: ");
    scanf("%f",&EnergyScale);

    printf("\nUse energy calibration? (1 = yes, 0 = no): ");
    scanf("%d",&ECal);

    char gainfile[150];

    if(ECal)
    {
        printf("\n Gain parameter file: ");
        scanf("%s", gainfile);
    }
    else
    {
        gainfile[0] = '\0';
    }

    int rootAnalysisWrite = 0;
    int rootTreeWrite = 0;
    int storeTrace = 0;

    printf("\nWrite ROOT Analysis (histograms, matrices, gates)? (1 = yes, 0 = no): ");
    scanf("%d", &rootAnalysisWrite);

    printf("\nWrite ROOT Tree (event data)? (1 = yes, 0 = no): ");
    scanf("%d", &rootTreeWrite);

    if(rootTreeWrite) {
        printf("\nStore trace data in ROOT Tree? (1 = yes, 0 = no): ");
        scanf("%d", &storeTrace);
    }

    // ========================================================
    // TIME SORTING TOGGLE
    // ========================================================

    int timeSortingChoice = 1;
    while(1)
    {
        printf("\nEnable TimeSorting? (1 = yes, 0 = no): ");
        if(scanf("%d", &timeSortingChoice) != 1)
        {
            printf("Invalid input. Enter 1 or 0.\n");
            while(getchar() != '\n');
            continue;
        }
        if(timeSortingChoice != 0 && timeSortingChoice != 1)
        {
            printf("Enter 1 to sort or 0 to use raw files.\n");
            continue;
        }
        break;
    }
    const bool TimeSorting = (timeSortingChoice == 1);

    // ========================================================
    // SELECT MODULE TYPE
    // ========================================================

    int modtype;

    while(1)
    {
        printf("Module type:\n");
        printf("1 = 100 MHz\n");
        printf("2 = 250 MHz\n");
        printf("3 = 500 MHz\n");
        printf("Choice: ");

        if(scanf("%d",&modtype) != 1)
        {
            printf("Invalid input. Enter integer.\n");

            while(getchar()!='\n');

            continue;
        }

        if(modtype < 1 || modtype > 3)
        {
            printf("Invalid choice.\n");
            continue;
        }

        break;
    }

    PixieType type;

    switch(modtype)
    {
        case 1:
            type = PixieType::PIXIE_100MHZ;
            break;

        case 2:
            type = PixieType::PIXIE_250MHZ;
            break;

        default:
            type = PixieType::PIXIE_500MHZ;
            break;
    }

    // ========================================================

    printf("Using %d OpenMP threads\n",
           omp_get_max_threads());

    // ========================================================
    // NUMBER OF FILES
    // ========================================================

    while(1)
    {
        printf("Number of files: ");

        if(scanf("%d",&Nfiles) != 1)
        {
            printf("Invalid input. Enter integer.\n");

            while(getchar()!='\n');

            continue;
        }

        if(Nfiles <= 0 || Nfiles > MAX_FILES)
        {
            printf("Enter value between 1 and %d\n",
                   MAX_FILES);

            continue;
        }

        break;
    }

    // ========================================================
    // FILE NAMES
    // ========================================================

    for(int i=0;i<Nfiles;i++)
    {
        while(1)
        {
            printf("File %d: ",i);

            if(scanf("%149s",file_list[i]) != 1)
            {
                printf("Invalid input.\n");

                while(getchar()!='\n');

                i--;

                continue;
            }

            FILE *test =
                fopen(file_list[i],"rb");

            if(!test)
            {
                printf("File does not exist.\n");
                continue;
            }

            fclose(test);

            break;
        }
    }

    // ========================================================
    // TIME WINDOW
    // ========================================================

    while(1)
    {
        printf("Time window (ns): ");

        if(scanf("%lf",&timewin) != 1)
        {
            printf("Invalid input. Enter number.\n");

            while(getchar()!='\n');

            continue;
        }

        if(timewin <= 0)
        {
            printf("Time window must be positive.\n");
            continue;
        }

        break;
    }

    // ========================================================

    // ========================================================
    // GROUP FILES BY RUN NUMBER
    // Filename convention expected by extract_run_number():
    // ..._<run-number>_mod<module>...
    // A run can contain any number of module files (up to MAX_FILES
    // input paths overall). All files with the same run number are
    // merged together; different runs are processed in parallel.
    // ========================================================

    std::map<int, std::vector<std::string>> grouped_runs;
    for(int i = 0; i < Nfiles; ++i)
    {
        const std::string name(file_list[i]);
        const int run = extract_run_number(name);
        if(run < 0)
        {
            std::cerr << "Cannot extract run number from filename: "
                      << name << "\n"
                      << "Expected a filename containing _<run-number>_mod.\n";
            return 1;
        }
        grouped_runs[run].push_back(name);
    }

    std::vector<int> run_numbers;
    std::vector<std::vector<std::string>> run_files;
    run_numbers.reserve(grouped_runs.size());
    run_files.reserve(grouped_runs.size());
    for(auto& item : grouped_runs)
    {
        run_numbers.push_back(item.first);
        run_files.push_back(std::move(item.second));
    }

    std::cout << "Grouped " << Nfiles << " input file(s) into "
              << run_numbers.size() << " run(s).\n";

    TimeScale = 1000.0 / TimeBinWidthPs;
    GainMatch(96, ECal, gainfile);

    Histograms globalH;

    allocate_histograms(&globalH);

    int nthreads = omp_get_max_threads();
    RootOutputManager rootOutput;
    if(!rootOutput.initialize(
            nthreads, type,
            rootAnalysisWrite != 0,
            rootTreeWrite != 0,
            storeTrace != 0)) {
        std::cerr << "Failed to initialize ROOT output.\n";
        free_histograms(&globalH);
        return 1;
    }

    // ========================================================
    // PARALLEL
    // ========================================================

    #pragma omp parallel
    {
        int tid =
            omp_get_thread_num();

        Histograms localH;

        allocate_histograms(&localH);

        // ====================================================
        // GROUP INPUT FILES BY RUN. Each run is one work item;
        // module files within a run are merged chronologically. TimeSorting
        // controls whether each module is externally sorted first or read raw.
        // ====================================================

        #pragma omp for schedule(dynamic, 1)
        for(int r = 0; r < (int)run_numbers.size(); ++r)
        {
            const int run_number = run_numbers[r];
            const std::vector<std::string>& modules = run_files[r];

            std::ostringstream msg;
            msg << "Thread " << tid << " processing run " << run_number
                << " using " << modules.size() << " module file(s):\n";
            for(const auto& name : modules)
                msg << "    " << name << "\n";

            #pragma omp critical
            std::cout << msg.str();

            unsigned int rng_state =
                static_cast<unsigned int>(time(nullptr) ^
                    (static_cast<unsigned int>(tid) * 0x9E3779B9u) ^
                    (static_cast<unsigned int>(run_number) * 0x85EBCA6Bu));

            process_file(
                modules,
                timewin,
                type,
                run_number,
                TimeSorting,
                &rng_state,
                [&](const CoincidenceBlock& block) {
                    // Both output managers consume the same precomputed addback
                    // values carried by CoincidenceBlock.
                    add_coincidence_to_histograms(block, &localH);

                    if(rootAnalysisWrite || rootTreeWrite) {
                        rootOutput.fillCoincidence(tid, run_number, block);

                        if(rootTreeWrite && storeTrace) {
                            for(const auto& member : block.members) {
                                if(!member.energyValid || member.event.trace.empty())
                                    continue;
                                const PixieEvent& ev = member.event;
                                rootOutput.fillTrace(tid, run_number, block.coincID,
                                    member.channel, ev.trace_length,
                                    member.calibratedEnergy, ev.coarse_time,
                                    ev.fine_time, ev.trace);
                            }
                        }
                    }
                }
            );
        }

        // ====================================================
        // MERGE HISTOGRAMS
        // ====================================================

        #pragma omp critical(merge_lock)
        {
            merge_histograms(
                &globalH,
                &localH
            );
        }

        free_histograms(&localH);
    }

    rootOutput.finalize();

    if(!write_binary_histograms(globalH)) {
        free_histograms(&globalH);
        return 1;
    }

    printf("Done.\n");
    free_histograms(&globalH);

    auto t_end =
    std::chrono::high_resolution_clock::now();

    double elapsed =
        std::chrono::duration<double>(
            t_end - t_start
        ).count();

    printf(
        "\nTotal elapsed time: %.3f seconds\n",
        elapsed
    );

    return 0;
}