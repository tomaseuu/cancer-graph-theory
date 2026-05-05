from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.api.routes import router
from app.db.database import DATABASE_ENABLED, init_db


app = FastAPI(
    title="OncoGraph API",
    description="Graph-based cancer gene discovery API",
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://127.0.0.1:3000",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.on_event("startup")
def startup_event():
    if DATABASE_ENABLED:
        init_db()


app.include_router(router)
