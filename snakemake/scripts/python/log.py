import datetime


class Log():
    def __init__(self):
        pass
    # 类方法
    @classmethod
    def getDatetime(cls):
        current_datetime = datetime.datetime.now()
        formatted_datetime = current_datetime.strftime("%m/%d/%Y %H:%M:%S")
        return formatted_datetime

    @classmethod
    def info(cls, text, file=None):
        if file:
            print(f"[{cls.getDatetime()}] INFO: {text}", file=file)
        else:
            print(f"[{cls.getDatetime()}] INFO: {text}")

    @classmethod
    def warning(cls, text, file):
        print(f"[{cls.getDatetime()}] WARNING: {text}", file=file)

