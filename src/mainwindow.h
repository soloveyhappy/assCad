#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QMdiArea>
#include <QDockWidget>
#include <QTreeView>

class QLabel;
class TreeModel;
class TreeItem;

class MainWindow : public QMainWindow
{
    Q_OBJECT
public:
    MainWindow();
    ~MainWindow();

protected:
    void mousePressEvent(QMouseEvent *) override;
signals:

public slots:
    void clicked(bool);
private slots:
    void onAction(bool);
private:
    QMdiArea* m_mdiArea;
    QDockWidget* m_dataDoc;
    QTreeView* m_treeView;
    TreeModel* treeModel;
    TreeItem* group1;
    TreeItem* group2;
    TreeItem* root;
};

#endif // MAINWINDOW_H
