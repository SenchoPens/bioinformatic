from typing import Callable, Any

def perform(condition: Callable[[], bool]):
    def wrap(action: Callable[[], Any]):
        if condition():
            return action()
        return None
    return wrap

def echo(condition: Callable[[], bool]):
    def wrap(msg: str, **kwargs):
        if condition():
            print(msg, **kwargs)
        return None
    return wrap
