#include <ncurses.h>
#include "3dcube.cpp"

int main() {
    initscr();
    start_color();
    curs_set(0);
    noecho();
    use_default_colors(); 

    int max_y, max_x;
    getmaxyx(stdscr, max_y, max_x);

    int d_min_max = min(max_x, max_y) / 2;
    int center_x = max_x / 2;
    int center_y = max_y / 2;

    std::vector<Vec> cubeVectors;
    generate_cube(cubeVectors, d_min_max);

    refresh();

    char gradient[] = " .:!/r(l1Z4H9W8@";
    double cube_max = 2.2;

    while (1) {
        clear();

        // Вращение куба
        rotateVectors(1, 1, 1, cubeVectors);
        sortVectors(cubeVectors);

        // Отрисовка куба
        for (const auto& vec : cubeVectors) {
            int col = 17 - static_cast<int>((vec.z + cube_max) / cube_max / 2 * 17);

            int screen_x = static_cast<int>((vec.x * d_min_max + center_x));
            int screen_y = static_cast<int>(vec.y * d_min_max * 11/24  + center_y);

            mvaddch(screen_y, screen_x, gradient[col]);
        }
        cout << endl;

        refresh();
        napms(25);
    }

    endwin();
    return 0;
}