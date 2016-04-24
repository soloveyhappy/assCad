#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include <QTextEdit>





class QLabel;

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    MainWindow();
    ~MainWindow();

protected:
    void mousePressEvent(QMouseEvent *);
signals:

public slots:
    void clicked(bool);
private:
    QDockWidget* dataDoc;
    QTextEdit* text;
};

#endif // MAINWINDOW_H
