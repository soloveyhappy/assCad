#include <QApplication>

#include "maincontroller.h"

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    int code = 0;
    {

        MainController controller;
        controller.show();
        code = app.exec();
    }
    return code;
}
