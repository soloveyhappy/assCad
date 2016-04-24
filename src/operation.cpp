#include "operation.h"

Operation::Operation(const QString& name)
    :m_name(name)
{

}

QString Operation::name() const
{
    return m_name;
}


Operation::~Operation()
{

}
