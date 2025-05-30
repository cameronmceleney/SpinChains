// The MIT License (MIT)
//
// Copyright (c) 2019 Luigi Pertoldi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// ============================================================================
//  ___   ___   ___   __    ___   ____  __   __   ___    __    ___
// | |_) | |_) / / \ / /`_ | |_) | |_  ( (` ( (` | |_)  / /\  | |_)
// |_|   |_| \ \_\_/ \_\_/ |_| \ |_|__ _)_) _)_) |_|_) /_/--\ |_| \_
//
// Very simple progress bar for c++ loops with internal running variable
//
// Author: Luigi Pertoldi
// Created: 3 dic 2016
//
// Notes: The bar must be used when there's no other possible source of output
//        inside the for loop
//

#ifndef __PROGRESSBAR_HPP
#define __PROGRESSBAR_HPP

#include <iostream>
#include <string>
#include <stdexcept>

class progressbar {

public:
    // default destructor
    ~progressbar()                             = default;

    // delete everything else
    progressbar           (progressbar const&) = delete;
    progressbar& operator=(progressbar const&) = delete;
    progressbar           (progressbar&&)      = delete;
    progressbar& operator=(progressbar&&)      = delete;

    // default constructor, must call set_niter later
    inline progressbar();
    inline explicit progressbar(int n, bool showbar=true);

    // reset bar to use it again
    inline void reset();
    // set number of loop iterations
    inline void set_niter(int iter);
    // chose your style
    inline void set_done_char(const std::string& sym) {done_char = sym;}
    inline void set_todo_char(const std::string& sym) {todo_char = sym;}
    inline void set_opening_bracket_char(const std::string& sym) {opening_bracket_char = sym;}
    inline void set_closing_bracket_char(const std::string& sym) {closing_bracket_char = sym;}
    // to show only the percentage
    inline void show_bar(bool flag = true) {do_show_bar = flag;}
    // main function
    inline void update();
    inline void updateIrregular();

private:
    int progress;
    int n_cycles;
    int last_perc;
    bool do_show_bar;
    bool update_is_called;

    std::string done_char;
    std::string todo_char;
    std::string opening_bracket_char;
    std::string closing_bracket_char;
};

inline progressbar::progressbar() :
        progress(0),
        n_cycles(0),
        last_perc(-1),
        do_show_bar(true),
        update_is_called(false),
        done_char("#"),
        todo_char(" "),
        opening_bracket_char("{"),
        closing_bracket_char("}") {}

inline progressbar::progressbar(int n, bool showbar) :
        progress(0),
        n_cycles(n),
        last_perc(0),
        do_show_bar(showbar),
        update_is_called(false),
        done_char("#"),
        todo_char(" "),
        opening_bracket_char("{"),
        closing_bracket_char("}") {}

inline void progressbar::reset() {
    progress = 0,
    update_is_called = false;
    last_perc = 0;
}

inline void progressbar::set_niter(int niter) {
    if (niter <= 0) throw std::invalid_argument(
                "progressbar::set_niter: number of iterations null or negative");
    n_cycles = niter + 1;
}

inline void progressbar::update() {

    if (n_cycles == 0) throw std::runtime_error(
                "progressbar::update: number of cycles not set");

    if (!update_is_called) {
        if (do_show_bar) {
            std::cout << opening_bracket_char;
            for (int _ = 0; _ < 50; _++) std::cout << todo_char;
            std::cout << closing_bracket_char << " 0%";
        }
        else std::cout << "0%";
    }
    update_is_called = true;

    int perc = 0;

    // compute percentage, if did not change, do nothing and return
    perc = static_cast<int>(progress*100.0/(n_cycles-1));
    if (perc < last_perc) return;

    // update percentage each unit
    if (perc == last_perc + 1) {
        // erase the correct  number of characters
        if      (perc <= 10)                std::cout << "\b\b"   << perc << '%';
        else if (perc  > 10 and perc < 100) std::cout << "\b\b\b" << perc << '%';
        else if (perc == 100)               std::cout << "\b\b\b" << perc << '%';
    }
    if (do_show_bar) {
        // update bar every ten units
        if (perc % 2 == 0) {
            // erase closing bracket
            std::cout << std::string(closing_bracket_char.size(), '\b');
            // erase trailing percentage characters
            if      (perc  < 10)               std::cout << "\b\b\b";
            else if (perc >= 10 && perc < 100) std::cout << "\b\b\b\b";
            else if (perc == 100)              std::cout << "\b\b\b\b\b";

            // erase 'todo_char'
            for (int j = 0; j < 50-(perc-1)/2; ++j) {
                std::cout << std::string(todo_char.size(), '\b');
            }

            // add one additional 'done_char'
            if (perc == 0) std::cout << todo_char;
            else           std::cout << done_char;

            // refill with 'todo_char'
            for (int j = 0; j < 50-(perc-1)/2-1; ++j) std::cout << todo_char;

            // read trailing percentage characters
            std::cout << closing_bracket_char << ' ' << perc << '%';
        }
    }
    last_perc = perc;
    ++progress;
    std::cout << std::flush;

    if (progress == 100) {
        // Small change I made to make completion more obvious to myself.
        std::cout << "\b\b\b\bFinished!";
    }
}

inline void progressbar::updateIrregular() {
        if (n_cycles == 0) throw std::runtime_error("progressbar::update: number of cycles not set");

        // compute percentage dynamically
        double perc = (static_cast<double>(progress) / n_cycles) * 100;

        // update progress bar and percentage display only if percentage changes
        if (static_cast<int>(perc) > last_perc) {
            std::cout << "\r" << opening_bracket_char;
            for (int i = 0; i < 50; ++i) {
                if (i < perc * 0.5) std::cout << done_char;
                else std::cout << todo_char;
            }
            std::cout << closing_bracket_char << ' ' << std::fixed << std::setprecision(1) << perc << '%';
            std::cout << std::defaultfloat;
            last_perc = static_cast<int>(perc);

            // check for completion
            if (progress == n_cycles) {
                std::cout << " Finished!";
            }

            std::cout << std::flush;
        }

        // increment progress
        ++progress;
    }

#endif