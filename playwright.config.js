const { defineConfig, devices } = require("@playwright/test");

const PORT = Number(process.env.BIOWORDLE_E2E_PORT || 4179);
const BASE_URL = `http://127.0.0.1:${PORT}`;

module.exports = defineConfig({
  testDir: "./tests/e2e",
  testMatch: /.*\.spec\.js$/,
  fullyParallel: false,
  workers: 1,
  retries: 0,
  timeout: 30_000,
  expect: {
    timeout: 8_000,
  },
  outputDir: "test-results",
  reporter: [
    ["list"],
    ["html", { open: "never", outputFolder: "playwright-report" }],
  ],
  use: {
    baseURL: BASE_URL,
    headless: true,
    trace: "retain-on-failure",
    screenshot: "only-on-failure",
    video: "off",
    viewport: { width: 1280, height: 900 },
  },
  projects: [
    {
      name: "chromium",
      use: { ...devices["Desktop Chrome"] },
    },
  ],
  webServer: {
    command: "node tests/e2e/static_server.mjs",
    url: `${BASE_URL}/index.html`,
    reuseExistingServer: false,
    timeout: 20_000,
    env: {
      ...process.env,
      BIOWORDLE_E2E_PORT: String(PORT),
    },
  },
});

