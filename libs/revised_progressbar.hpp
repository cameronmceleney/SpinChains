#ifndef __PROGRESSBAR_HPP
#define __PROGRESSBAR_HPP

#include <iostream>
#include <string>
#include <iomanip> // For std::setw and std::setfill
#include <stdexcept>

class progressbar {

public:
    progressbar(int n = 100, bool showbar = true);
    ~progressbar() = default;

    void reset();
    void set_niter(int iter);
    void set_done_char(const std::string& sym) { done_char = sym; }
    void set_todo_char(const std::string& sym) { todo_char = sym; }
    void set_opening_bracket_char(const std::string& sym) { opening_bracket_char = sym; }
    void set_closing_bracket_char(const std::string& sym) { closing_bracket_char = sym; }
    void show_bar(bool flag = true) { do_show_bar = flag; }
    void update();
    void update_with_iterations(int current, int total);
    void set_progress_flags(bool* globalFlag, bool* localFlag);

private:
    int progress, n_cycles, last_perc;
    bool do_show_bar, update_is_called;
    bool *isUsingProgressBar, *hasLocalProgressBar; // Pointers to flags
    std::string done_char, todo_char, opening_bracket_char, closing_bracket_char;

    void print_progress(int perc, int current = -1, int total = -1);
};

progressbar::progressbar(int n, bool showbar) : progress(0), n_cycles(n), last_perc(0), do_show_bar(showbar), update_is_called(false), isUsingProgressBar(nullptr), hasLocalProgressBar(nullptr), done_char("#"), todo_char(" "), opening_bracket_char("{"), closing_bracket_char("}") {}

void progressbar::reset() {
    progress = 0;
    last_perc = 0;
    update_is_called = false;
}

void progressbar::set_niter(int niter) {
    if (niter <= 0) throw std::invalid_argument("progressbar::set_niter: number of iterations null or negative");
    n_cycles = niter;
}

void progressbar::update() {
    ++progress;
    int perc = static_cast<int>(100.0 * progress / n_cycles);
    if (perc > last_perc) {
        print_progress(perc);
        last_perc = perc;
    }
}

void progressbar::update_with_iterations(int current, int total) {
    int perc = static_cast<int>(100.0 * current / total);
    if (perc > last_perc || current == total) { // Ensure final update
        print_progress(perc, current, total);
        last_perc = perc;
    }
}

void progressbar::set_progress_flags(bool* globalFlag, bool* localFlag) {
    isUsingProgressBar = globalFlag;
    hasLocalProgressBar = localFlag;
}

void progressbar::print_progress(int perc, int current, int total) {
    if (do_show_bar) {
        std::cout << '\r' << opening_bracket_char; // Use '\r' to return cursor to the start of the line
        int pos = 50 * perc / 100;
        for (int i = 0; i < 50; ++i) {
            if (i < pos) std::cout << done_char;
            else std::cout << todo_char;
        }
        std::cout << closing_bracket_char << ' ';
    }
    else std::cout << '\r'; // Ensure we're starting from the beginning of the line

    if (current >= 0 && total > 0) { // Handling iteration tracking
        std::cout << std::setw(3) << perc << "% " << current << "/" << total;
    }
    else { // Default percentage-only mode
        std::cout << perc << "%";
    }
    std::cout << std::flush; // Ensure output is immediately visible
}

#endif // __PROGRESSBAR_HPP