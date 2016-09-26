#include <QApplication>
#include <QDebug>
#include "maincontroller.h"

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    QTranslator translator;
    translator.load("qt_" + QLocale::system().name(), QLibraryInfo::location(QLibraryInfo::TranslationsPath));
    app.setTranslator(&translator);
    int code = 0;
    {

        MainController controller;
        controller.show();
        code = app.exec();
    }
    return code;
}
